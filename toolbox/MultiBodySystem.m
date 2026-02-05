classdef MultiBodySystem  < handle
    properties
        Name (1,1) string = "UnnamedSystem"
        gravity (3,1) sym = 0
        bodies struct = struct()                % All bodies in the system
        time (1,1) sym = sym('time', 'real')
        dof struct = struct()                   % Degrees of freedom
        dof_idx struct = struct()        
        q (:,1) sym = []
        doc struct = struct()                   % Degrees of constraint
        doc_idx struct = struct()               % Constraintforces are calculated in the direction of these coordinates
        z (:,1) sym = []                        % Constraint coordinates
        inputs struct = struct()
        externals struct = struct()             % Externally calculated value that may depend on inputs and states
        external_deps struct = struct()         % Dependencies of the external
        external_params (:,1) string = []

        outputs struct = struct()

        params Parameters = []

        attached_bodies (1,:) Body = Body.empty % Bodies attached to the inertial frame

        aux_state struct = struct()             % auxiliary state names
        aux_impl_ode (:,1) sym = []             % auxiliary implicit first order ode        

        setupCompleted = false                  % set to true once model is finished an no further data can be added

        eom (:,1) sym = []
        Fz (:,1) sym = []                       % cash for calculated constraint forces
        M (:,:) sym = []
        C (:,:) sym = []
        K (:,:) sym = []
        B (:,:) sym = []
        CD (:,:) sym = []
        F (:,:) sym = []

        keep_dof logical = []                   % cache for positional states that are not unused
        aux_ode_first_order logical = []        % cache for test result, if all aux ODEs are first order
    end

    properties (Access = private)
    end

    methods
        % Constructor with optional name
        function obj = MultiBodySystem(name, dof, input)
            obj.params = Parameters();
            if nargin > 0
                obj.Name = string(name);
            end
            if nargin>1
                obj.addGeneralizedCoordinate(dof)
            end
            if nargin>2
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
            obj.checkName(coordName)
    
            % Define symbolic variable dynamically
            save_assumptions = assumptions;
            symFun = str2sym([coordName '(time)']);
            assumptions(save_assumptions); % str2sym deletes old assumtions
            assumeAlso(symFun, 'real')
    
            % Store it
            obj.dof.(coordName) = symFun;
            obj.dof.([coordName '_d']) = diff(symFun, obj.time);
            obj.dof.([coordName '_dd']) = diff(symFun, obj.time, 2);
            obj.q(end+1,1) = symFun;
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
            obj.checkName(coordName)
    
            % Define symbolic variable dynamically
            save_assumptions = assumptions;
            symFun = str2sym([coordName '(time)']);
            assumptions(save_assumptions); % str2sym deletes the original assumptions
            assumeAlso(symFun, 'real')

            % Store it
            obj.doc.(coordName) = symFun;
            obj.doc.([coordName '_d']) = diff(symFun, obj.time);
            obj.doc.([coordName '_dd']) = diff(symFun, obj.time, 2);
            obj.z(end+1,1) = symFun;
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
    
            % Define symbolic variable dynamically
            symVar = sym(inName, 'real');
            assume(symVar, 'real')
    
            % Store it
            obj.inputs.(inName) = symVar;
        end

        % Add input by name
        function addOutput(obj, outName, expr)
            arguments
                obj
                outName { mustBeTextScalar }
                expr (1,1) sym
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
                depends (1,:) sym = []
            end
    
            obj.checkSetupNotCompleted();            
            checkName(obj, extName)
    
            % Define symbolic variable dynamically
            if isempty(depends)
                symFun = sym(extName, 'real');
            else
                save_assumptions = assumptions;
                symFun = str2sym(string(extName) + "(" + join(string(depends), ',') + ")");
                assumptions(save_assumptions); % str2sym deletes the original assumptions
                assumeAlso(symFun, 'real')
            end
    
            % Store it
            obj.externals.(extName) = symFun;
            obj.external_deps.(extName) = depends;
        end

        % Add auxilliary state by name
        function addAuxState(obj, auxName)
            arguments
                obj
                auxName { MultiBodySystem.mustBeNonemptyCharOrCell }
            end
            if iscell(auxName)
                for i = 1:length(auxName)
                    obj.addAuxState(auxName{i})
                end
                return
            end
    
            obj.checkSetupNotCompleted();
            obj.checkName(auxName)
    
            % Define symbolic variable dynamically
            save_assumptions = assumptions;
            symFun = str2sym([auxName '(time)']);
            assumptions(save_assumptions); % str2sym deletes the original assumptions
            assumeAlso(symFun, 'real')

            % Store it
            obj.aux_state.(auxName) = symFun;
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
                var (:,1) sym
                order (1,1) double = 1
            end

            dx = diff(var, obj.time, order);
        end

        function x = getAuxState(obj, ii)
            fn = fieldnames(obj.aux_state);
            x = sym.empty;
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
            obj.checkAuxODEFirstOrder()

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

        function n = getQName(obj, i, naming)
            arguments
                obj 
                i (1,:) double = 1:length(obj.q)
                naming {mustBeText} = 'real_name'
            end
            if isempty(i)
                i = 1:length(obj.q);
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
            if nargin>1 && isscalar(i)
                n = n{1};
            end
        end

        function n = getQdName(obj, i, naming)
            arguments
                obj 
                i (1,:) double = 1:length(obj.q)
                naming {mustBeText} = 'real_name'
            end
            if isempty(i)
                i = 1:length(obj.q);
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
            if isscalar(n)
                n = n{1};
            end
        end

        function n = getQddName(obj, i, naming)
            arguments
                obj 
                i (1,:) double = 1:length(obj.q)
                naming {mustBeText} = 'real_name'
            end
            if isempty(i)
                i = 1:length(obj.q);
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
            if isscalar(n)
                n = n{1};
            end
        end

        function n = getAuxName(obj, i, deriv, naming)
            arguments
                obj 
                i (1,:) double = []
                deriv (1,1) double = 0
                naming {mustBeText} = 'real_name'
            end
            if isempty(i)
                i = 1:obj.getNumAux();
            end
            if deriv>0
                deriv_str = repmat('d', 1, deriv);
            else
                deriv_str = 0;
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
            if isscalar(n)
                n = n{1};
            end
        end

        function n = getStateNames(obj, eliminate_unused)
            arguments
                obj 
                eliminate_unused (1,1) logical = false
            end
            
            obj.checkAuxODEFirstOrder()
            n = obj.getQName();
            if eliminate_unused
                keep_dof_ = obj.getUsedPositionalStates();
                n = n(keep_dof_);
            end
            
            n_qd = obj.getQdName();
            if ~iscell(n_qd)
                n_qd = {n_qd};
            end
            n_aux = obj.getAuxName();
            if ~iscell(n_aux)
                n_aux = {n_aux};
            end
            n = [n; n_qd; n_aux];
        end

        function nd = getDStateNames(obj, eliminate_unused)
            arguments
                obj 
                eliminate_unused (1,1) logical = false
            end
            
            obj.checkAuxODEFirstOrder()
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
                expr sym
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

        function [partial_names, deps, dep_names] = getExternalDerivs(obj, i, naming)
            fn = fieldnames(obj.externals);
            n = fn{i};
            deps = obj.external_deps.(n);
            dep_names = obj.replaceDOFs(deps);
            partial_names = cell(1, length(deps));
            for i = 1:length(deps)
                partial_names{i} = ['d' n '_d' char(dep_names(i))];
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
                eom_ = sym(zeros(length(obj.q), 1));
                for i= 1:length(obj.attached_bodies)
                    eom_ = eom_ + obj.attached_bodies(i).collectGenForces;
                end
    
                eom_ = simplify(eom_, 'Steps', 50);
                eom_ = collect(eom_, struct2array(obj.dof));
                obj.eom = eom_;
            end
            if with_aux
                eom_ = [eom_; obj.aux_impl_ode];
            end
        end

        function [f_impl, i_state_idx] = getImplStateSpaceODE(obj, eliminate_unused)
            arguments
                obj 
                eliminate_unused (1,1) logical = false
            end
            
            % For now, only first order auxilliary ODEs are permitted
            obj.checkAuxODEFirstOrder()

            if eliminate_unused
                keep_dof_ = obj.getUsedPositionalStates();
            else
                keep_dof_ = true(length(obj.q), 1);
            end

            for i = 1:length(obj.q)
                % positonal state derivatives
                % these have to match with the names produced in getDStateNames
                dx1(i, 1) = sym(sprintf('dot_%s', obj.getQName(i)), 'real');
            end

            eom_ = getEOM(obj);

            % TODO: check if some auxilliary ODEs are also intgrals
            i_state_idx = [true(sum(keep_dof_), 1); false(length(eom_) + length(obj.aux_impl_ode), 1)];

            f_impl = [dx1(keep_dof_) - obj.getTimeDeriv(obj.q(keep_dof_), 1) ; eom_; obj.aux_impl_ode];
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
                Fz_ = sym(zeros(length(obj.z), 1));
                for i= 1:length(obj.attached_bodies)
                    Fz_ = Fz_ + obj.attached_bodies(i).collectConstrForces;
                end
                
                Fz_ = simplify(Fz_);
                obj.Fz = Fz_;
            end

            Fz_ = Fz_(ismember(fieldnames(obj.doc_idx), name));
            if remove_eps
                Fz_ = Body.removeEps(Fz_);
            end
        end

        function fun = eomFunO2(obj, filename)
            eom_ = getEOM(obj);
            [eom_, vars]= obj.replaceVars(eom_, 'numbered');
            if ~exist('filename', 'var')
                fun = matlabFunction(eom_, 'Vars', {vars.q_, vars.qd_, vars.qdd_, struct2array(obj.inputs), struct2array(obj.params.getParamSyms())});            
            else
                matlabFunction(eom_, 'File', filename, 'Vars', {vars.q_, vars.qd_, vars.qdd_, struct2array(obj.inputs), struct2array(obj.params.getParamSyms())});
            end
        end

        function fun = eomFunDescriptor(obj, filename)
            eom_ =  getEOM(obj);
            [eom_, vars]= obj.replaceVars(eom_, 'numbered');
            M_ = -getM(obj);
            f_= simplify(eom_ + M_*vars.qdd);

            if ~exist('filename', 'var')
                fun = matlabFunction(M_, f_, 'Vars', {vars.q_, vars.qd_, struct2array(obj.inputs), struct2array(obj.params.getParamSyms())});            
            else
                matlabFunction(M_, f_, 'File', filename, 'Vars', {vars.q_, vars.qd_, struct2array(obj.inputs), struct2array(obj.params.getParamSyms())});
            end
        end

        function fun = eomDae(obj, filename)
            eom_ = getEOM(obj);
            
            if isempty(fieldnames(obj.inputs))
                inputs_ = sym('dummy', 'real');
            else
                inputs_ = struct2array(obj.inputs);
            end
            % replace external functions by simple symbols because
            % daeFunction doesn't work with symfuns
            fn = fieldnames(obj.externals);
            if isempty(fn)
                externalSyms = sym('dummy', 'real');
            else
                externalSyms= sym.empty(1, 0);
                for i = 1:length(fn)
                    externalSyms(end+1) = sym(fn{i}, 'real');
                end
                eom_ = subs(eom_, struct2array(obj.externals), externalSyms);
            end

            for i = 1:length(obj.q)
                save_assumptions = assumptions;
                x1(i, 1) = symfun(str2sym(sprintf('x1_%d(time)', i)), obj.time);
                x2(i, 1) = symfun(str2sym(sprintf('x2_%d(time)', i)), obj.time);
                assumptions(save_assumptions); % str2sym deletes the original assumptions
                assumeAlso(x1(i, 1), 'real')
                assumeAlso(x2(i, 1), 'real')
            end
            x = [x1 ; x2; struct2array(obj.aux_state)];

            eom_ = subs(eom_, diff(obj.q, obj.time, 2), diff(x2, obj.time));
            eom_ = subs(eom_, diff(obj.q, obj.time), x2);
            eom_ = subs(eom_, obj.q, x1);

            f_impl = [diff(x1, obj.time) == x2 ; eom_; obj.aux_impl_ode];

            if ~exist('filename', 'var')
                fun = daeFunction(f_impl, x, inputs_, externalSyms, struct2array(obj.params.getParamSyms()));
            else
                daeFunction(f_impl, x, inputs_, externalSyms, struct2array(obj.params.getParamSyms()), 'File', filename);
            end
        end

        function [e, vars] = replaceDOFs(obj, e, naming)
            arguments
                obj 
                e (:,:) sym
                naming {mustBeText} = 'real_name'
            end
            % TODO: save assumptions here too?
            vars.qdd = str2sym(obj.getQddName([], naming));
            assume(vars.qdd, 'real')
            vars.qd = str2sym(obj.getQdName([], naming));
            assume(vars.qd, 'real')
            vars.q = str2sym(obj.getQName([], naming));
            assume(vars.q, 'real')

            e = subs(e, diff(obj.q, obj.time, 2), vars.qdd);
            e = subs(e, diff(obj.q, obj.time), vars.qd);
            e = subs(e, obj.q, vars.q);

            vars.aux = str2sym(obj.getAuxName([], 0, naming));
            assume(vars.aux, 'real')
            vars.auxd = str2sym(obj.getAuxName([], 1, naming));
            assume(vars.aux, 'real')
            vars.auxdd = str2sym(obj.getAuxName([], 2, naming));
            assume(vars.aux, 'real')

            e = subs(e, obj.getTimeDeriv(struct2array(obj.aux_state), 2), vars.auxdd);
            e = subs(e, obj.getTimeDeriv(struct2array(obj.aux_state), 1), vars.auxd);
            e = subs(e, struct2array(obj.aux_state)', vars.aux);
        end

        function [e, vars] = replaceVars(obj, e, naming)
            arguments
                obj 
                e (:,:) sym
                naming {mustBeText} = 'real_name'
            end
            % external derivatives have to be replaced first to not get confused when
            % the DOFs and externals are replaced
            exts = struct2array(obj.externals);
            ext_names = fieldnames(obj.externals);
            ext_d = cell(size(exts));
            for i = 1:length(exts)
                partial_names = obj.getExternalDerivs(i, naming);
                % TODO: save assumptions here too?
                partial_syms = str2sym(partial_names);
                assume(partial_syms, 'real')
                deps = obj.external_deps.(ext_names{i});
                for j = 1:length(deps)
                    % partial_deriv = functionalDerivative(obj.externals(i), deps(j));
                    partial_deriv = diff(exts(i), deps(j));
                    e = subs(e, partial_deriv, partial_syms(j));            
                end
                ext_d{i} = partial_syms;
            end
            ext_vars = str2sym(obj.getExternalName([], naming));
            assume(ext_vars, 'real')
            e = subs(e, exts, ext_vars');

            [e, vars] = replaceDOFs(obj, e, naming);
            
            vars.ext = ext_vars;
            vars.ext_d = ext_d;

            % TODO: save_assumptions here too?
            vars.u = str2sym(obj.getInName([], naming));
            assume(vars.u, 'real')
            vars.p = str2sym(obj.getParamName([], naming));
            assume(vars.p, 'real')


            e = subs(e, struct2array(obj.inputs), vars.u');
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
                vars (:, 1) sym
                cached sym
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
            if with_aux
                vars = [vars; obj.getTimeDeriv(struct2array(obj.aux_state), 2)];
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
                    vars = [obj.q; struct2array(obj.aux_state); obj.getTimeDeriv(obj.q, 1); obj.getTimeDeriv(struct2array(obj.aux_state), 1); struct2array(obj.inputs)'];
                else
                    vars = [obj.q; obj.getTimeDeriv(obj.q, 1); struct2array(obj.inputs)'];
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
                    vars = [vars; obj.getTimeDeriv(struct2array(obj.aux_state), 2)];
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

        function name_in_use= checkName(obj, name)
            if nargout>0
                name_in_use= false;
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
            if isempty(obj.keep_dof)
                eom_ = obj.getEOM();
                % we need to remove the diff so as not to find the
                % functions inside them
                dummy = @(y) sym('dummy');
                eom_ = mapSymType(eom_, 'diff', dummy);
                outs= struct2array(obj.outputs);
                if isempty(outs)
                    time_funs = findSymType(eom_, 'symfun');
                else
                    outs = mapSymType(outs, 'diff', dummy);
                    time_funs = [findSymType(eom_, 'symfun') findSymType(outs, 'symfun')];
                end

                keep_dof_ = true(length(obj.q), 1);
                for i = 1:length(obj.q)
                    if ~ismember(obj.q(i), time_funs)
                        keep_dof_(i) = false;
                    end
                end
                obj.keep_dof = keep_dof_;
            else
                keep_dof_ = obj.keep_dof;
            end
        end

        function checkAuxODEFirstOrder(obj)
            obj.checkSetupCompleted()
            if isempty(obj.aux_ode_first_order)
                % get all derivatives in the aux ODEs
                derivs_aux = findSymType(obj.aux_impl_ode, 'diff');
                % for each check if it is of an aux state and if yes the order
                for i = 1:length(derivs_aux)
                    c = children(derivs_aux(i));
                    % derivatimes should be only of individual time dependent
                    % functions, these are then the first child of diff
                    if ismember(c{1}, struct2array(obj.aux_state))
                        % for now assume, all diffs are of time only, thus the
                        % number of children determines the order
                        if length(c)>2
                            obj.aux_ode_first_order = false;
                            error('One or more auxilliary ODEs is not of first order.')
                        end
                    end
                end
                obj.aux_ode_first_order = true;
            else
                if ~obj.aux_ode_first_order
                    error('One or more auxilliary ODEs is not of first order.')
                end
            end
        end
    end
    
    methods (Static, Access = private)
        function mustBeNonemptyCharOrCell(v)
            if ~(ischar(v) && ~isempty(v) || iscellstr(v) && ~isempty(v))
                error('Input must be a nonempty char vector or nonempty cell array of char vectors.');
            end
        end
    end
end
