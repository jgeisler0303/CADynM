classdef MultiBodySystem  < handle
    properties
        Name (1,1) string = "UnnamedSystem"
        gravity
        bodies struct = struct()                % All bodies in the system
        time
        dof struct = struct()                   % Degrees of freedom
        dof_idx struct = struct()        
        q
        doc struct = struct()                   % Degrees of constraint
        doc_idx struct = struct()               % Constraintforces are calculated in the direction of these coordinates
        z                                       % Constraint coordinates
        inputs struct = struct()
        externals struct = struct()             % Externally calculated value that may depend on inputs and states
        external_deps struct = struct()         % Dependencies of the external
        external_params (:,1) string = []

        outputs struct = struct()

        params Parameters

        attached_bodies (1,:) Body = Body.empty % Bodies attached to the inertial frame

        aux_state struct = struct()             % auxiliary state names
        aux_order struct = struct()             % order = max derivative of each state, currently only first order is supported
        aux_impl_ode                            % auxiliary implicit first order ode        

        setupCompleted = false                  % set to true once model is finished an no further data can be added

        eom
        Fz                                      % cache for calculated constraint forces
        M
        C 
        K 
        B 
        CD 
        F

        keep_positional_states logical = []     % cache for positional states that are not unused
        aux_ode_first_order logical = []        % cache for test result, if all aux ODEs are first order

        sym_eps                                 % variable to control cross terms of small (elastic) deformations
        sym_eps_rot                             % variable to control cross terms of small (elastic) rotations
        dummy
    end

    methods
        % Constructor with optional name and symbolic backend
        % symbolicBackend: 'sym' for MATLAB Symbolic Toolbox, 'msym' for MAMaS/Maxima (default)
        function obj = MultiBodySystem(name, dof, input)
            % Initialize symbolic variables based on selected backend
            % we have to assign it here, otherwise it may be the same for all instances
            obj.gravity = obj.sym([0 0 0]');
            obj.time = obj.sym('time');
            obj.sym_eps = obj.sym('eps');
            obj.sym_eps_rot = obj.sym('eps_rot');
            obj.dummy = obj.sym('dummy');  

            obj.q = obj.sym([]);
            obj.z = obj.sym([]);
            obj.aux_impl_ode = obj.sym([]);

            % Pass system reference to Parameters so it can use the correct backend
            obj.params = Parameters([], obj);
            if nargin > 0 && ~isempty(name)
                obj.Name = string(name);
            end
            if nargin > 1 && ~isempty(dof)
                obj.addGeneralizedCoordinate(dof)
            end
            if nargin > 2 && ~isempty(input)
                obj.addInput(input)
            end
        end

        % Add a generalized coordinate by name
        function addGeneralizedCoordinate(obj, coordName)
            arguments
                obj
                coordName { MultiBodySystem.mustBeNonemptyCharOrCell }
            end

            if iscell(coordName)
                for i = 1:length(coordName)
                    obj.addGeneralizedCoordinate(coordName{i})
                end
                return
            end

            obj.checkSetupNotCompleted();

            [obj.dof.(coordName), obj.dof.([coordName '_d']), obj.dof.([coordName '_dd'])] = createDerivatives(obj, coordName);
    
            obj.q(end+1,1) = obj.dof.(coordName);
            obj.dof_idx.(coordName) = length(obj.q);
        end

        % Add a generalized constraint coordinate by name
        function addConstraintCoordinate(obj, coordName)
            arguments
                obj
                coordName { MultiBodySystem.mustBeNonemptyCharOrCell }
            end

            obj.checkSetupNotCompleted();

            if iscell(coordName)
                for i = 1:length(coordName)
                    obj.addConstraintCoordinate(coordName{i})
                end
                return
            end
    
            obj.checkSetupNotCompleted();
            [obj.doc.(coordName), obj.doc.([coordName '_d']), obj.doc.([coordName '_dd'])] = createDerivatives(obj, coordName);
    
            obj.z(end+1,1) = obj.doc.(coordName);
            obj.doc_idx.(coordName) = length(obj.z);
        end

        % Add a parameter by name
        function p = addParameter(obj, paramName, dims, value)
            arguments
                obj
                paramName
                dims = []
                value = []
            end
            obj.checkSetupNotCompleted();

            if isstruct(paramName)
                obj.params.setParamRefStruct(paramName)
            elseif iscell(paramName)
                p_ = {};
                for i = 1:length(paramName)
                    if ~iscell(dims)                
                        dims_ = dims;
                    else
                        dims_ = dims{i};
                    end
                    if ~iscell(value)                
                        value_ = value;
                    else
                        value_ = value{i};
                    end
                    p_{end+1} = obj.addParameter(paramName{i}, dims_, value_);
                end
            else
                p_ = obj.params.addParameter(paramName, dims, value);
            end

            if nargout>0
                p = p_;
            end
        end

        function p = addExternalParameter(obj, paramName, dims, value)
            arguments
                obj
                paramName
                dims = []
                value = []
            end
            p = addParameter(obj, paramName, dims, value);
            if isstruct(paramName)
                paramName = fieldnames(paramName);
            end
            obj.external_params = [obj.external_params; string(paramName)];
        end

        function setParamValue(obj, name, value)
            obj.params.addParameter(name, [], value)
        end

        % Add input by name
        function addInput(obj, inName)
            arguments
                obj
                inName { MultiBodySystem.mustBeNonemptyCharOrCell }
            end

            if iscell(inName)
                for i = 1:length(inName)
                    obj.addInput(inName{i})
                end
                return
            end
    
            obj.checkSetupNotCompleted();
            obj.checkName(inName)
    
            % Define symbolic variable dynamically using backend abstraction
            symVar = obj.sym(inName);
            % TODO: consider time derivatives?
    
            % Store it
            obj.inputs.(inName) = symVar;
        end

        % Add output by name
        function addOutput(obj, outName, expr)
            arguments
                obj
                outName { mustBeTextScalar }
                expr
            end
            
            if isfield(obj.outputs, outName)
                error('Output "%s" already defined.', outName)
            end
            obj.outputs.(outName) = expr;
        end

        % Add external value by name
        function addExternal(obj, extName, depends)
            arguments
                obj
                extName (1,:) char
                depends = []
            end
    
            obj.checkSetupNotCompleted();            
            checkName(obj, extName)
    
            % Define symbolic variable dynamically using backend abstraction
            symFun = obj.sym(extName, depends);
    
            % Store it
            obj.externals.(extName) = symFun;
            obj.external_deps.(extName) = depends;
        end

        % Add auxilliary state by name
        function addAuxState(obj, auxName, order)
            arguments
                obj
                auxName { MultiBodySystem.mustBeNonemptyCharOrCell }
                order (:,1) double = 1
            end
            if iscell(auxName)
                for i = 1:length(auxName)
                    obj.addAuxState(auxName{i}, order(i))
                end
                return
            end
            if order>1
                error('Only first order auxilliary ODEs are currently supported.')
            end

            obj.checkSetupNotCompleted();
            obj.checkName(auxName)
    
            % Define symbolic variable dynamically using backend abstraction
            % TODO: use getDerivatives here too, somehow manage number of used derivatives
            symFun = obj.sym(auxName, obj.time);

            % Store it
            obj.aux_state.(auxName) = symFun;
            obj.aux_order.(auxName) = order;
        end

        % Add auxilliary state by name
        function addAuxImplODE(obj, ode)
            obj.checkSetupNotCompleted();
            obj.aux_impl_ode(end+1) = ode;
        end

        % Add a body with a unique name
        function addBody(obj, body)
            arguments
                obj
                body (1,1) Body
            end
    
            obj.checkSetupNotCompleted();
            if isempty(body.Name)
                name_in_use= true;
                i= 1;
                while name_in_use
                    new_name= sprintf('body%02d', i);
                    name_in_use= checkName(new_name);
                    i = i+1;
                end
                body.Name = new_name;
            else
                obj.checkName(body.Name)
            end

            obj.bodies.(body.Name) = body;
        end

        % Add a direct child
        function addChild(obj, body)
            arguments
                obj
                body (1,1) Body
            end

            obj.checkSetupNotCompleted();
            obj.addBody(body);
            
            body.parent= obj;
            body.system= obj;
            obj.attached_bodies(end+1)= body;
        end

        function completeSetup(obj)
            obj.checkSetupNotCompleted()
            obj.setupCompleted = true;
            
            obj.prepareKinematics();
        end

        function removeUnusedParameters(obj)
            vars = [symvar(obj.eom) symvar(struct2array(obj.outputs)) symvar(obj.aux_impl_ode)];
            obj.params.removeUnused([string(vars) obj.external_params(:)']);
        end

        % Get number of bodies
        function n = getNumBodies(obj)
            n = numel(fieldnames(obj.bodies));
        end

        % Display all bodies
        function displaySystem(obj)
            fprintf('Multi-Body System: %s\n', obj.Name);
            fn = fieldnames(obj.bodies);
            fprintf('Number of bodies: %d\n', numel(fn));
            for i = 1:numel(fn)
                fprintf('  %d. %s (%s)\n', i, fn{i}, class(obj.bodies.(fn{i})));
            end
        end
        
        % Get the time derivative of auxilliary state
        function dx = getTimeDeriv(obj, var, order)
            arguments
                obj
                var
                order (1,1) double = 1
            end
            aux_var_idx = ismember(struct2array(obj.aux_state), var);
            if any(aux_var_idx)
                aux_order_vec = struct2array(obj.aux_order);
                if any(aux_order_vec(aux_order_vec)<order)
                    error('Requested derivative order exceeds defined order for this auxilliary state for (one of) "%s".', string(var))
                end
            end
            dx = diff(var, obj.time, order);
        end

        function x = getAuxState(obj, ii)
            fn = fieldnames(obj.aux_state);
            x = obj.sym([]);
            for i = ii
                x(end+1) = obj.aux_state.(fn{i});
            end
        end

        function n = getNumDOF(obj)
            n = length(obj.q);
        end
        
        % TODO also add getNumAuxStates for higher order auxilliary ODEs
        function n = getNumAux(obj)
            n = length(fieldnames(obj.aux_state));
        end

        function n = getNumStates(obj, eliminate_unused)
            arguments
                obj 
                eliminate_unused (1,1) logical = false
            end

            if eliminate_unused
                keep_dof_ = obj.getUsedPositionalStates();
            else
                keep_dof_ = true(length(obj.q), 1);
            end
            
            n = sum(keep_dof_) + obj.getNumDOF + obj.getNumAux();
        end

        function n = getNumIn(obj)
            n = length(fieldnames(obj.inputs));
        end

        function n = getNumOut(obj)
            n = length(fieldnames(obj.outputs));
        end

        function n = getNumParams(obj)
            n = obj.params.getNumParams();
        end

        function n = getNumExternals(obj)
            n = length(fieldnames(obj.externals));
        end

        function n = getQName(obj, i_, naming)
            arguments
                obj 
                i_ (1,:) double = []
                naming {mustBeText} = 'real_name'
            end
            if isempty(i_)
                i = 1:length(obj.q);
            else
                i = i_;
            end
            switch naming
                case 'real_name'
                    fn = fieldnames(obj.dof_idx);
                    n = fn(i);
                case 'numbered'
                    n = sprintfc('q_%d', i)';
                case 'cpp'
                    n = sprintfc('q_IDX%dXDI_', i-1)';
            end
            if nargin>1 && isscalar(i_)
                n = n{1};
            end
        end

        function n = getQdName(obj, i_, naming)
            arguments
                obj 
                i_ (1,:) double = []
                naming {mustBeText} = 'real_name'
            end
            if isempty(i_)
                i = 1:length(obj.q);
            else
                i = i_;
            end
            switch naming
                case 'real_name'
                    fn = fieldnames(obj.dof_idx);
                    n = strcat(fn(i), '_d');
                case 'numbered'
                    n = sprintfc('qd_%d', i)';
                case 'cpp'
                    n = sprintfc('qd_IDX%dXDI_', i-1)';
            end
            if isscalar(n) && isscalar(i_)
                n = n{1};
            end
        end

        function n = getQddName(obj, i_, naming)
            arguments
                obj 
                i_ (1,:) double = 1:length(obj.q)
                naming {mustBeText} = 'real_name'
            end
            if isempty(i_)
                i = 1:length(obj.q);
            else
                i = i_;
            end
            switch naming
                case 'real_name'
                    fn = fieldnames(obj.dof_idx);
                    n = strcat(fn(i), '_dd');
                case 'numbered'
                    n = sprintfc('qdd_%d', i)';
                case 'cpp'
                    n = sprintfc('qdd_IDX%dXDI_', i-1)';
            end
            if isscalar(n) && isscalar(i_)
                n = n{1};
            end
        end

        function n = getAuxName(obj, i_, deriv, naming)
            arguments
                obj 
                i_ (1,:) double = []
                deriv (1,1) double = 0
                naming {mustBeText} = 'real_name'
            end
            if isempty(i_)
                i = 1:obj.getNumAux();
            else
                i = i_;
            end
            if deriv>0
                if deriv>2
                    deriv_str = ['d' num2str(deriv)];
                else
                    deriv_str = repmat('d', 1, deriv);
                end
            else
                deriv_str = '';
            end

            switch naming
                case 'real_name'
                    fn = fieldnames(obj.aux_state);
                    if deriv>0
                        n = strcat(fn(i), ['_' deriv_str]);
                    else
                        n = fn(i);
                    end
                case 'numbered'
                    basename = 'aux';
                    if deriv>0
                        basename = [basename deriv_str];
                    end
                    n = sprintfc('%s_%d', basename, i)';
                case 'cpp'
                    n = sprintfc(['q' deriv_str '_IDX%dXDI_'], i-1+obj.getNumDOF)';
            end
            if isscalar(n) && isscalar(i_)
                n = n{1};
            end
        end

        function n = getStateNames(obj, eliminate_unused)
            arguments
                obj 
                eliminate_unused (1,1) logical = false
            end
            
            n = obj.getQName();
            if eliminate_unused
                keep_dof_ = obj.getUsedPositionalStates();
                n = n(keep_dof_);
            end
            
            n_qd = obj.getQdName();
            n_aux = obj.getAuxName();
            n = [n; n_qd; n_aux];
        end

        function nd = getDStateNames(obj, eliminate_unused)
            arguments
                obj 
                eliminate_unused (1,1) logical = false
            end
            
            nd = obj.getQName();
            if eliminate_unused
                keep_dof_ = obj.getUsedPositionalStates();
                nd = nd(keep_dof_);
            end
            nd = strcat('dot_', nd);

            n_qdd = obj.getQddName();
            if ~iscell(n_qdd)
                n_qdd = {n_qdd};
            end
            n_auxd = obj.getAuxName([], 1);
            if ~iscell(n_auxd)
                n_auxd = {n_auxd};
            end
            nd = [nd; n_qdd; n_auxd];
        end

        function n = getInName(obj, i, naming)
            arguments
                obj 
                i (1,:) double = []
                naming {mustBeText} = 'real_name'
            end
            if isempty(i)
                i = 1:length(fieldnames(obj.inputs));
            end
            switch naming
                case 'real_name'
                    fn = fieldnames(obj.inputs);
                    n = fn(i);
                case 'numbered'
                    n = sprintfc('in_%d', i)';
                case 'cpp'
                    n = sprintfc('u_IDX%dXDI_', i-1)';
            end
            if isscalar(n)
                n = n{1};
            end
        end

        function n = getOutName(obj, i)
            fn = fieldnames(obj.outputs);
            n = fn{i};
        end

        function n = getParamName(obj, i, naming)
            arguments
                obj 
                i (1,:) double = []
                naming {mustBeText} = 'real_name'
            end
            if isempty(i)
                i = 1:obj.params.getNumParams();
            end
            switch naming
                case 'real_name'
                    fn = obj.params.getParamNames();
                    n = fn(i)';
                case 'numbered'
                    n = sprintfc('p_%d', i)';
                case 'cpp'
                    fn = obj.params.getParamNames();
                    n = strcat('PSTRUCT_', fn(i))';
            end
            if isscalar(n)
                n = n{1};
            end
        end

        function n = getExternalName(obj, i, naming)
            arguments
                obj 
                i (1,:) double = []
                naming {mustBeText} = 'real_name'
            end
            if isempty(i)
                i = 1:length(fieldnames(obj.externals));
            end
            switch naming
                case 'real_name'
                    fn = fieldnames(obj.externals);
                    n = fn(i);
                case 'numbered'
                    n = sprintfc('ext_%d', i)';
                case 'cpp'
                    fn = fieldnames(obj.externals);
                    n = fn(i);
            end
            if isscalar(n)
                n = n{1};
            end
        end

        function tf = isConstant(obj, expr, naming)
            arguments
                obj
                expr
                naming {mustBeTextScalar} = 'real_name'
            end
        
            constants_list = string(obj.getParamName([], naming));
            syms_in_expr = symvar(expr);
        
            % If no symbolic variables, it is constant
            if isempty(syms_in_expr)
                tf = true;
                return;
            end
        
            % Expression is constant if all detected symbols are allowed
            tf = all(ismember(string(syms_in_expr), constants_list));
        end

        function [partial_names, deps, dep_names] = getExternalDerivs(obj, i, ~)
            fn = fieldnames(obj.externals);
            n = fn{i};
            deps = obj.external_deps.(n);
            dep_names = obj.replaceDOFs(deps);
            partial_names = cell(1, length(deps));
            for j = 1:length(deps)
                partial_names{j} = ['d' n '_d' char(dep_names(j))];
            end
        end

        function eom_ = getEOM(obj, with_aux)
            arguments
                obj 
                with_aux (1,1) logical = false
            end
            obj.checkSetupCompleted()
            if ~isempty(obj.eom)
                eom_ = obj.eom;
            else
                eom_ = obj.sym(zeros(length(obj.q), 1));
                for i= 1:length(obj.attached_bodies)
                    eom_ = eom_ + obj.attached_bodies(i).collectGenForces;
                end
    
                eom_ = obj.simplify(eom_);
                obj.eom = eom_;
            end
            if with_aux && ~isempty(obj.aux_impl_ode)
                % work-aroud for inconsistent MAMaS concatenation rules: pre allocate correct dimensions
                eom__ = obj.sym(zeros(length(eom_)+length(obj.aux_impl_ode), 1));
                eom__(1:length(eom_), 1) = eom_;
                eom__(length(eom_)+1:end, 1) = obj.aux_impl_ode;
                eom_ = eom__;
            end
        end

        function [f_impl, i_state_idx] = getImplStateSpaceODE(obj, eliminate_unused)
            arguments
                obj 
                eliminate_unused (1,1) logical = false
            end
            
            if eliminate_unused
                keep_dof_ = obj.getUsedPositionalStates();
            else
                keep_dof_ = true(length(obj.q), 1);
            end

            dx1 = obj.sym([]);
            for i = 1:length(obj.q)
                % positonal state derivatives
                % these have to match with the names produced in getDStateNames
                dx1(i, 1) = obj.sym(sprintf('dot_%s', obj.getQName(i)));
            end

            eom_ = getEOM(obj);

            % TODO: check if some auxilliary ODEs are also intgrals
            i_state_idx = [true(sum(keep_dof_), 1); false(length(eom_) + length(obj.aux_impl_ode), 1)];

            % work-aroud for inconsistent MAMaS concatenation rules: pre allocate correct dimensions
            f_impl = obj.sym(zeros(sum(keep_dof_)+obj.getNumDOF+length(obj.aux_impl_ode), 1));

            f_impl(1:sum(keep_dof_)) = dx1(keep_dof_) - obj.getTimeDeriv(obj.q(keep_dof_), 1);
            f_impl((1:obj.getNumDOF)+sum(keep_dof_)) = eom_;
            f_impl((1:length(obj.aux_impl_ode))+sum(keep_dof_)+obj.getNumDOF) = obj.aux_impl_ode;
        end

        function Fz_ = getConstraintForce(obj, name, remove_eps)
            arguments
                obj 
                name 
                remove_eps = true;
            end
            obj.checkSetupCompleted()
            if ~isempty(obj.eom)
                Fz_ = obj.Fz;
            else
                Fz_ = obj.sym(zeros(length(obj.z), 1));
                for i= 1:length(obj.attached_bodies)
                    Fz_ = Fz_ + obj.attached_bodies(i).collectConstrForces;
                end
                
                Fz_ = obj.simplify(Fz_);
                obj.Fz = Fz_;
            end

            Fz_ = Fz_(ismember(fieldnames(obj.doc_idx), name));
            if remove_eps
                Fz_ = obj.removeEps(Fz_);
            end
        end

        function [e, vars] = replaceDOFs(obj, e, naming)
            arguments
                obj 
                e (:,:)
                naming {mustBeText} = 'real_name'
            end
            % TODO: save assumptions here too?
            vars.qdd = cellfun(@(n)obj.sym(n), obj.getQddName([], naming));
            vars.qd = cellfun(@(n)obj.sym(n), obj.getQdName([], naming));
            vars.q = cellfun(@(n)obj.sym(n), obj.getQName([], naming));

            e = subs(e, diff(obj.q, obj.time, 2), vars.qdd);
            e = subs(e, diff(obj.q, obj.time), vars.qd);
            e = subs(e, obj.q, vars.q);

            vars.aux = cellfun(@(n)obj.sym(n), obj.getAuxName([], 0, naming));
            vars.auxd = cellfun(@(n)obj.sym(n), obj.getAuxName([], 1, naming));
            vars.auxdd = cellfun(@(n)obj.sym(n), obj.getAuxName([], 2, naming));
            
            % TODO: handle arbitrary order aux states
            e = subs(e, obj.getTimeDeriv(struct2array(obj.aux_state), 1), vars.auxd);
            e = subs(e, struct2array(obj.aux_state).', vars.aux);
        end

        function [e, vars] = replaceVars(obj, e, naming)
            arguments
                obj 
                e (:,:) 
                naming {mustBeText} = 'real_name'
            end
            % external derivatives have to be replaced first to not get confused when
            % the DOFs and externals are replaced
            exts = struct2array(obj.externals);
            ext_names = fieldnames(obj.externals);
            ext_d = cell(size(exts));
            for i = 1:length(exts)
                partial_names = obj.getExternalDerivs(i, naming);
                partial_syms = cellfun(@(n)obj.sym(n), partial_names);
                deps = obj.external_deps.(ext_names{i});
                for j = 1:length(deps)
                    partial_deriv = diff(exts(i), deps(j));
                    e = subs(e, partial_deriv, partial_syms(j));            
                end
                ext_d{i} = partial_syms;
            end
            ext_vars = cellfun(@(n)obj.sym(n), obj.getExternalName([], naming));
            e = subs(e, exts, ext_vars.');

            [e, vars] = replaceDOFs(obj, e, naming);
            
            vars.ext = ext_vars;
            vars.ext_d = ext_d;

            vars.u = cellfun(@(n)obj.sym(n), obj.getInName([], naming));
            vars.p = cellfun(@(n)obj.sym(n), obj.getParamName([], naming));

            e = subs(e, struct2array(obj.inputs), vars.u.');
            e = subs(e, struct2array(obj.params.getParamSyms()), vars.p);
        end

        function ap = getParameterArray(obj)
            fn = obj.params.getParamNames();
            ap = [];
            for i = 1:length(fn)
                if ~obj.params.hasValue(fn{i})
                    error('No value set for parameter "%s".', fn{i})
                end
                p = obj.params.getValue(fn{i});
                ap = [ap p(:)];
            end
        end

        function J = getJacobian(obj, vars, cached, with_aux)
            arguments
                obj 
                vars (:, 1) 
                cached 
                with_aux (1,1) logical = false
            end
            obj.checkSetupCompleted()
            n_states = obj.getNumDOF;
            if with_aux
                n_states = n_states + obj.getNumAux;
            end
            if size(cached, 1)>=n_states
                J = cached(1:n_states, 1:n_states);
            else
                eom_ = obj.getEOM(with_aux);

                J= jacobian(eom_, vars);
            end
        end
        % TODO: refactor to use one common jacobian function called with
        % different derivative vectors
        % Generalized mass matrix
        function M_ = getM(obj, with_aux)
            arguments
                obj 
                with_aux (1,1) logical = false
            end
            vars = obj.getTimeDeriv(obj.q, 2);
            % currently only first order aux odes are allowed
            if with_aux
                aux_names = fieldnames(obj.aux_order);
                aux_state_dd = cell(length(aux_names), 1);
                for i = 1:length(aux_names)
                    if obj.aux_order.(aux_names{i})<2
                        % for second order integrators like Newmark Beta we
                        % need to hanlde each ode like a second order ode
                        aux_state_dd{i} = obj.dummy;
                    elseif obj.aux_order.(aux_names{i})<3
                        aux_state_dd{i} = obj.getTimeDeriv(obj.aux_state.(aux_names{i}), 2);
                    end
                end
                vars = [vars; cell2mat(aux_state_dd)];
            end

            M_ = getJacobian(obj, vars, obj.M, with_aux);
            obj.M = M_;
        end

        % Generalized coriolis and damping matrix
        function C_ = getC(obj, with_aux)
            arguments
                obj 
                with_aux (1,1) logical = false
            end
            vars = obj.getTimeDeriv(obj.q, 1);
            if with_aux
                vars = [vars; obj.getTimeDeriv(struct2array(obj.aux_state), 1)];
            end

            C_ = getJacobian(obj, vars, obj.C, with_aux);
            obj.C = C_;
        end

        % Generalized stiffness matrix
        function K_ = getK(obj,  with_aux)
            arguments
                obj 
                with_aux (1,1) logical = false
            end
            vars = obj.q;
            if with_aux
                vars = [vars; struct2array(obj.aux_state)];
            end

            K_ = getJacobian(obj, vars, obj.K, with_aux);
            obj.K = K_;
        end

        % Input Jacobian
        function B_ = getB(obj, with_aux)
            arguments
                obj 
                with_aux (1,1) logical = false
            end
            obj.checkSetupCompleted()
            n_eqns = obj.getNumDOF;
            if with_aux
                n_eqns = n_eqns + obj.getNumAux;
            end
            if isempty(fieldnames(obj.inputs)) || size(obj.B, 1)>=n_eqns
                B_ = obj.B(1:n_eqns, :);
            else
                vars = struct2array(obj.inputs);
                B_= jacobian(obj.getEOM(with_aux), vars);
                obj.B = B_;
            end
        end

        % Output Jacobian
        function CD_ = getCD(obj, with_aux)
            arguments
                obj 
                with_aux (1,1) logical = false
            end
            obj.checkSetupCompleted()

            n_derivs = 2*obj.getNumDOF + obj.getNumIn;
            if with_aux
                n_derivs = n_derivs + 2*obj.getNumAux;
            end
            % CD may hold the jacobian with or without auxilliaries
            % if it was prior computed without, it has to be recomputed
            if isempty(fieldnames(obj.outputs)) || size(obj.CD, 2)>=n_derivs
                if isempty(fieldnames(obj.outputs)) || size(obj.CD, 2)==n_derivs
                    % it was last computed with auxilliaries and is
                    % demanded again or it was computed without and is
                    % demanded as such: no problem
                    CD_ = obj.CD;
                else
                    % it was last computed with and is no demanded without
                    % auxilliaries
                    idx = [true(1, obj.getNumDOF) false(1, obj.getNumAux)];
                    idx = [idx idx true(1, obj.getNumIn)];
                    CD_ = obj.CD(:, idx);
                end
            else
                if with_aux
                    vars = [obj.q; struct2array(obj.aux_state); obj.getTimeDeriv(obj.q, 1); obj.getTimeDeriv(struct2array(obj.aux_state), 1); struct2array(obj.inputs).'];
                else
                    vars = [obj.q; obj.getTimeDeriv(obj.q, 1); struct2array(obj.inputs).'];
                end
                CD_= jacobian(struct2array(obj.outputs), vars);
                obj.CD = CD_;
            end
        end

        % Output Jacobian wrt accelerations
        function F_ = getF(obj, with_aux)
            arguments
                obj 
                with_aux (1,1) logical = false
            end
            obj.checkSetupCompleted()
            n_derivs = obj.getNumDOF;
            if with_aux
                n_derivs = n_derivs + obj.getNumAux;
            end
            if isempty(fieldnames(obj.outputs)) || size(obj.F, 2)>=n_derivs
                F_ = obj.F(:, 1:n_derivs);
            else
                vars = obj.getTimeDeriv(obj.q, 2);
                if with_aux
                    aux_names = fieldnames(obj.aux_order);
                    aux_state_dd = cell(length(aux_names), 1);
                    for i = 1:length(aux_names)
                        if obj.aux_order.(aux_names{i})<2
                            % for second order integrators like Newmark Beta we
                            % need to hanlde each ode like a second order ode
                            aux_state_dd{i} = obj.dummy;
                        elseif obj.aux_order.(aux_names{i})<3
                            aux_state_dd{i} = obj.getTimeDeriv(obj.aux_state.(aux_names{i}), 2);
                        end
                    end
                    vars = [vars; cell2mat(aux_state_dd)];
                end
                F_= jacobian(struct2array(obj.outputs), vars);
                obj.F = F_;
            end
        end

        function v = paramVec(obj, p)
            pv = struct2array(obj.params.getParamSyms());
            v = zeros(size(pv));

            for i= 1:length(pv)
                v(i) = p.(char(pv(i)));
            end
        end

        function expr = removeEps(obj, expr, keep_symbols)
            arguments
                obj
                expr
                keep_symbols (1,1) logical = false
            end
            if obj.isMSym
                % remove small rotations
                expr = expr.subs(obj.sym_eps_rot^2, 0, true);
                if ~keep_symbols
                    expr = expr.subs(obj.sym_eps_rot, 1);
                end

                % remove cross terms of small elastic terms
                expr = expr.subs(obj.sym_eps^2, 0, true);
                if ~keep_symbols
                    expr = expr.subs(obj.sym_eps, 1);
                end
            elseif obj.isSym
                % remove small rotations
                expr = obj.removeHigherOrderTermsSym(expr, obj.sym_eps_rot, keep_symbols);

                % remove cross terms of small elastic terms
                % TODO: add reference why only first order terms are kept
                expr = obj.removeHigherOrderTermsSym(expr, obj.sym_eps, keep_symbols);
            else
                error('Unsupported symbolic backend: %s', obj.getsetSymbolicBackend)
            end
            % TODO: what about eps*eps_rot cross terms?
        end

        function expr = keepEps(obj, expr)
            expr = subs(expr, [obj.sym_eps, obj.sym_eps_rot], [1, 1]);
        end

        function expr = removeDOC(obj, expr)
            if obj.isMSym   
                expr = subs(expr, struct2array(obj.doc), 0);
            elseif obj.isSym
                expr = subs(expr, obj.z, zeros(size(obj.z)));
            else
                error('Unsupported symbolic backend: %s', obj.getsetSymbolicBackend)
            end
            for i = 1:length(obj.z)
                expr = subs(expr, obj.z(i), 0);
            end
        end

        function checkSetupCompleted(obj)
            if ~obj.setupCompleted
                error('Model setup has not been completed. Please run completeSetup first.')
            end
        end

        function checkSetupNotCompleted(obj)
            if obj.setupCompleted
                error('Model setup has been completed. No further data can be added.')

            end
        end
    end

    methods (Access = private)
        % calculate kinematics
        function prepareKinematics(obj)
            for i= 1:length(obj.attached_bodies)
                obj.attached_bodies(i).T0= obj.attached_bodies(i).T;
                obj.attached_bodies(i).prepareKinematics;
            end
        end
        
        function [x, x_d, x_dd] = createDerivatives(obj, var_name)
            var_name_d = [var_name '_d']; % This should be the sam as 'real_name' like in getQdName
            var_name_dd = [var_name '_dd'];

            obj.checkName(var_name);
            obj.checkName(var_name_d);
            obj.checkName(var_name_dd);
            if obj.isMSym
                x = obj.sym(var_name);
                gradef(x, obj.time, var_name_d)
                x_d  = obj.sym(var_name_d);
                gradef(x_d, obj.time, var_name_dd)
                x_dd = obj.sym(var_name_dd, obj.time);
            elseif obj.isSym
                x = str2sym([var_name '(' char(obj.time) ')']);
                x_d = diff(x, obj.time);
                x_dd = diff(x, obj.time, 2);
            else
                error('Unsupported symbolic backend: %s', obj.getsetSymbolicBackend)
            end
        end

        function name_in_use= checkName(obj, name)
            if nargout>0
                name_in_use= false;
            end
            if ismember(name, {'sym_eps', 'sym_eps_rot', 'dummy'})
                if nargout>0
                    name_in_use= true;
                else
                    error("The names sym_eps, sym_eps_rot and dummy are reserved for internal use.");
                end                
            end
            if isfield(obj.bodies, name)
                if nargout>0
                    name_in_use= true;
                else
                    error("Body name '%s' already exists in the system.", name);
                end
            end
            if isfield(obj.dof, name)
                if nargout>0
                    name_in_use= true;
                else
                    error("Coordinate name '%s' already exists in the system.", name);
                end
            end
            if isfield(obj.doc, name)
                if nargout>0
                    name_in_use= true;
                else
                    error("Constraint coordinate name '%s' already exists in the system.", name);
                end
            end
            % if obj.params.paramInUse(name)
            %     if nargout>0
            %         name_in_use= true;
            %     else
            %         error("Parameter name '%s' already exists in the system.", name);
            %     end
            % end
            if isfield(obj.inputs, name)
                if nargout>0
                    name_in_use= true;
                else
                    error("Input name '%s' already exists in the system.", name);
                end
            end
            if isfield(obj.externals, name)
                if nargout>0
                    name_in_use= true;
                else
                    error("External value name '%s' already exists in the system.", name);
                end
            end
            if isfield(obj.aux_state, name)
                if nargout>0
                    name_in_use= true;
                else
                    error("Auxilliary state name '%s' already exists in the system.", name);
                end
            end
        end

        function keep_dof_ = getUsedPositionalStates(obj)
            if isempty(obj.keep_positional_states)
                eom_ = obj.getEOM();
                outs= struct2array(obj.outputs);

                keep_dof_ = true(length(obj.q), 1);
                if obj.isSym
                    % we need to remove the diff so as not to find the
                    % functions inside them
                    dummy_fun = @(y) obj.dummy;
                    eom_ = mapSymType(eom_, 'diff', dummy_fun);
                    if isempty(outs)
                        time_funs = findSymType(eom_, 'symfun');
                    else
                        outs = mapSymType(outs, 'diff', dummy_fun);
                        time_funs = [findSymType(eom_, 'symfun') findSymType(outs, 'symfun')];
                    end

                    for i = 1:length(obj.q)
                        if ~ismember(obj.q(i), time_funs)
                            keep_dof_(i) = false;
                        end
                    end
                elseif obj.isMSym
                    for i = 1:length(obj.q)
                        if ~eom_.contains(obj.q(i)) && (isempty(outs) || ~outs.contains(obj.q(i)))
                            keep_dof_(i) = false;
                        end
                    end
                else
                    error('Unknown symbolic backend: %s', getSymbolicBackend());
                end

                obj.keep_positional_states = keep_dof_;
            else
                keep_dof_ = obj.keep_positional_states;
            end
        end
    end

    methods
        % ===== SYMBOLIC BACKEND ABSTRACTION METHODS =====
        % These methods provide a unified interface to work with symbolic objects
        % that work with both 'sym' (MATLAB Symbolic) and 'msym' (MAMaS/Maxima) backends.
        
        % Create a symbolic object based on the selected backend
        function result = sym(obj, val, deps)
            % sym(name, ...) or sym(value)
            % Delegates to the appropriate backend (sym or msym)
            if obj.isSym
                if nargin<3 || isempty(deps)
                    result = sym(val);
                else
                    result = str2sym(string(val) + "(" + join(string(deps), ',') + ")");
                end
            elseif obj.isMSym
                if nargin<3 || isempty(deps)
                    if val == pi
                        result = msym.pi;
                    else
                        result = msym(val);
                    end
                else
                    result = msym(val);
                    result.depends(deps);
                    % TODO: use gradef here too?
                end
            else
                error('Unknown symbolic backend: %s', getSymbolicBackend());
            end
        end

        function result = simplify(obj, expr)
            if obj.isSym
                result = simplify(expr, 'Steps', 50);
                result = collect(result, struct2array(obj.dof));
            elseif obj.isMSym
                result = expr.simplify();
            else
                error('Unknown symbolic backend: %s', getSymbolicBackend());
            end
        end

        function [opt, tnames, tvalues] = optimize(obj, varargin)
            if obj.isSym
                [tvalues, opt_, tnames] = feval2sym(symengine, 'symobj::optimizeWithIntermediates', [varargin{:}]);
                opt_ = opt_(1);

                col1 = 0;
                opt = cell(1, length(varargin));
                for i = 1:length(varargin)
                    opt{i} = opt_(:, (1:size(varargin{i}, 2))+col1);
                    col1 = col1 + size(varargin{i}, 2);
                end
            elseif obj.isMSym
                for i = 1:length(varargin)
                    if ~isscalar_matlab(varargin{i})
                        varargin{i} = msym.matrix(varargin{i});
                    end
                end
                [opt, tnames, tvalues] = optimize(cell2mat(varargin));
                if isscalar(varargin)
                    opt = {opt};
                else
                    opt = num2cell(opt);
                end
            else
                error('Unknown symbolic backend: %s', getSymbolicBackend());
            end
        end
    end


    methods (Static)
        function b = getsetSymbolicBackend(b)
            persistent symbolicBackend
            if nargin<1 || isempty(b)
                if isempty(symbolicBackend)
                    symbolicBackend = 'sym';
                end
            else
                symbolicBackend = b;
            end
            b = symbolicBackend;
        end

        function setSym()
            MultiBodySystem.getsetSymbolicBackend('sym');
        end
        
        function setMSym()
            MultiBodySystem.getsetSymbolicBackend('msym');
        end

        function tf = isSym()
            tf = strcmp(MultiBodySystem.getsetSymbolicBackend(), 'sym');
        end

        function tf = isMSym()
            tf = strcmp(MultiBodySystem.getsetSymbolicBackend(), 'msym');
        end
    end

    methods (Static, Access = private)
        function expr = removeHigherOrderTermsSym(expr, monom, keep_symbols)
            for i = 1:length(expr)
                c = coeffs(expr(i), monom, 'All');
                switch length(c)
                    case {0,1} % empty or no terms with monom, or only first order terms: keep them
                    otherwise % remove all but first order terms
                        if keep_symbols
                            expr(i) = c(end-1)*monom + c(end);
                        else
                            expr(i) = sum(c(end-1:end)); % eliminate monom
                        end
                end
            end
        end

        function mustBeNonemptyCharOrCell(v)
            if ~(ischar(v) && ~isempty(v) || iscellstr(v) && ~isempty(v))
                error('Input must be a nonempty char vector or nonempty cell array of char vectors.');
            end
        end
    end
end
