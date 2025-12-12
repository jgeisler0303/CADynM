classdef MultiBodySystem  < handle
    properties
        Name (1,1) string = "UnnamedSystem"
        gravity (3,1) sym = 0
        bodies struct = struct()  % Stores bodies by name
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

        outputs struct = struct()

        params struct = struct()

        children (1,:) Body = Body.empty

        aux_state struct = struct()             % auxiliary state names
        aux_impl_ode (:,1) sym = []             % auxiliary implicit first order ode        

        setupCompleted = false                  % set to true once model is finished an no further data can be added

        eom (:,1) sym = []
        M (:,:) sym = []
        C (:,:) sym = []
        K (:,:) sym = []
        B (:,:) sym = []
        CD (:,:) sym = []
        F (:,:) sym = []
    end

    properties (Access = private)
    end

    methods
        % Constructor with optional name
        function obj = MultiBodySystem(name, dof, input)
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

            obj.checkSetupCompleted();
            obj.checkName(coordName)
    
            % Define symbolic variable dynamically
            % symFun = symfun(str2sym([coordName '(time)']), obj.time);
            symFun = str2sym([coordName '(time)']);
    
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

            obj.checkSetupCompleted();

            if iscell(coordName)
                for i = 1:length(coordName)
                    obj.addConstraintCoordinate(coordName{i})
                end
                return
            end
    
            obj.checkSetupCompleted();
            obj.checkName(coordName)
    
            % Define symbolic variable dynamically
            % symFun = symfun(str2sym([coordName '(time)']), obj.time);
            symFun = str2sym([coordName '(time)']);
    
            % Store it
            obj.doc.(coordName) = symFun;
            obj.doc.([coordName '_d']) = diff(symFun, obj.time);
            obj.doc.([coordName '_dd']) = diff(symFun, obj.time, 2);
            obj.z(end+1,1) = symFun;
            obj.doc_idx.(coordName) = length(obj.z);
        end

        % Add a parameter by name
        function addParameter(obj, paramName, dims)
            arguments
                obj
                paramName { MultiBodySystem.mustBeNonemptyCharOrCell }
                dims (1,:) double = []
            end
            if iscell(paramName)
                if ~isempty(dims)
                    error('Non scalar Parameters cannot be created at once.')
                end
                for i = 1:length(paramName)
                    obj.addParameter(paramName{i})
                end
                return
            end
    
            obj.checkSetupCompleted();
            obj.checkName(paramName)
    
            % Define symbolic variable dynamically
            if ~isempty(dims)
                symVar = sym(paramName, dims, 'real');
            else
                symVar = sym(paramName, 'real');
            end
    
            % Store it
            obj.params.(paramName) = symVar;
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
    
            obj.checkSetupCompleted();
            obj.checkName(inName)
    
            % Define symbolic variable dynamically
            symVar = sym(inName, 'real');
    
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

            obj.checkSetupCompleted();
            obj.checkName(outName)

            obj.outputs.(outName) = expr;
        end

        % Add external value by name
        function addExternal(obj, extName, depends)
            arguments
                obj
                extName (1,:) char
                depends sym
            end
    
            obj.checkSetupCompleted();            
            checkName(obj, extName)
    
            % Define symbolic variable dynamically
            symFun = str2sym(string(extName) + "(" + join(string(depends), ',') + ")");
    
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
    
            obj.checkSetupCompleted();
            obj.checkName(auxName)
    
            % Define symbolic variable dynamically
            symFun = str2sym([auxName '(time)']);
    
            % Store it
            obj.aux_state.(auxName) = symFun;
        end

        % Add auxilliary state by name
        function addAuxImplODE(obj, ode)
            obj.checkSetupCompleted();
            obj.aux_impl_ode(end+1) = ode;
        end

        % Add a body with a unique name
        function addBody(obj, body)
            arguments
                obj
                body (1,1) Body
            end
    
            obj.checkSetupCompleted();
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

            obj.checkSetupCompleted();
            obj.addBody(body);
            
            body.parent= obj;
            body.system= obj;
            obj.children(end+1)= body;
        end

        function completeSetup(obj)
            obj.checkSetupCompleted()
            obj.setupCompleted = true;
            
            obj.prepareKinematics();
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
        function dx = getTimeDeriv(obj, dof, order)
            arguments
                obj
                dof (1,1) sym
                order (1,1) double = 1
            end

            dx = diff(dof, obj.time, order);
        end

        function n = getNumDOF(obj)
            n = length(obj.q);
        end

        function n = getNumStates(obj)
            % TODO: later change to possibly number reduced by unnecessary states
            n = 2*length(obj.q);
        end

        function n = getNumIn(obj)
            n = length(fieldnames(obj.inputs));
        end

        function n = getNumOut(obj)
            n = length(fieldnames(obj.outputs));
        end

        function n = getNumParams(obj)
            n = length(fieldnames(obj.params));
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
            if isscalar(n)
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
                i = 1:length(fieldnames(obj.params));
            end
            switch naming
                case 'real_name'
                    fn = fieldnames(obj.params);
                    n = fn(i)';
                case 'numbered'
                    n = sprintfc('p_%d', i)';
                case 'cpp'
                    fn = fieldnames(obj.params);
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
        %ISPURECONSTANT Determine whether a symbolic expression is constant.
        %   TF = ISPURECONSTANT(EXPR, CONSTANTS_LIST) returns TRUE if EXPR contains
        %   only:
        %       - numeric literals,
        %       - symbols in CONSTANTS_LIST,
        %   regardless of what symbolic functions are applied to them.
        %
        %   Thus, sin(a), log(a + pi), exp(3*C1) are constant if 'a' and 'C1'
        %   are included in CONSTANTS_LIST.
        
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

        function eom_ = getEOM(obj)
            if ~isempty(obj.eom)
                eom_ = obj.eom;
                return
            end
            eom_ = sym(zeros(length(obj.q), 1));
            for i= 1:length(obj.children)
                eom_ = eom_ + obj.children(i).collectGenForces;
            end
            
            eom_ = simplify(eom_);
            obj.eom = eom_;
        end

        function Fz = getConstraintForces(obj)
            Fz = sym(zeros(length(obj.z), 1));
            for i= 1:length(obj.children)
                Fz = Fz + obj.children(i).collectConstrForces;
            end
            
            Fz = simplify(Fz);
        end

        function fun = eomFunO2(obj, filename)
            eom_ = getEOM(obj);
            [eom_, vars]= obj.replaceVars(eom_, 'numbered');
            if ~exist('filename', 'var')
                fun = matlabFunction(eom_, 'Vars', {vars.q_, vars.qd_, vars.qdd_, struct2array(obj.inputs), struct2array(obj.params)});            
            else
                matlabFunction(eom_, 'File', filename, 'Vars', {vars.q_, vars.qd_, vars.qdd_, struct2array(obj.inputs), struct2array(obj.params)});
            end
        end

        function fun = eomFunDescriptor(obj, filename)
            eom_ =  getEOM(obj);
            [eom_, vars]= obj.replaceVars(eom_, 'numbered');
            M_ = -getM(obj, eom_);
            f_= simplify(eom_ + M_*vars.qdd_);

            if ~exist('filename', 'var')
                fun = matlabFunction(M_, f_, 'Vars', {vars.q_, vars.qd_, struct2array(obj.inputs), struct2array(obj.params)});            
            else
                matlabFunction(M_, f_, 'File', filename, 'Vars', {vars.q_, vars.qd_, struct2array(obj.inputs), struct2array(obj.params)});
            end
        end

        function fun = eomDae(obj, filename)
            eom_ = getEOM(obj);
            
            if isempty(fieldnames(obj.inputs))
                inputs_ = sym('dummy');
            else
                inputs_ = struct2array(obj.inputs);
            end
            % replace external functions by simple symbols because
            % daeFunction doesn't work with symfuns
            fn = fieldnames(obj.externals);
            if isempty(fn)
                externalSyms = sym('dummy');
            else
                externalSyms= sym.empty(1, 0);
                for i = 1:length(fn)
                    externalSyms(end+1) = sym(fn{i});
                end
                eom_ = subs(eom_, struct2array(obj.externals), externalSyms);
            end

            for i = 1:length(obj.q)
                x1(i, 1) = symfun(str2sym(sprintf('x1_%d(time)', i)), obj.time);
                x2(i, 1) = symfun(str2sym(sprintf('x2_%d(time)', i)), obj.time);
            end
            x = [x1 ; x2; struct2array(obj.aux_state)];

            eom_ = subs(eom_, diff(obj.q, obj.time, 2), diff(x2, obj.time));
            eom_ = subs(eom_, diff(obj.q, obj.time), x2);
            eom_ = subs(eom_, obj.q, x1);

            f_impl = [diff(x1, obj.time) == x2 ; eom_; obj.aux_impl_ode];

            if ~exist('filename', 'var')
                fun = daeFunction(f_impl, x, inputs_, externalSyms, struct2array(obj.params));
            else
                daeFunction(f_impl, x, inputs_, externalSyms, struct2array(obj.params), 'File', filename);
            end
        end

        function [e, vars] = replaceDOFs(obj, e, naming)
            arguments
                obj 
                e (:,:) sym
                naming {mustBeText} = 'real_name'
            end
            vars.qdd = str2sym(obj.getQddName([], naming));
            vars.qd = str2sym(obj.getQdName([], naming));
            vars.q = str2sym(obj.getQName([], naming));

            e = subs(e, diff(obj.q, obj.time, 2), vars.qdd);
            e = subs(e, diff(obj.q, obj.time), vars.qd);
            e = subs(e, obj.q, vars.q);
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
                [partial_names, deps] = obj.getExternalDerivs(i, naming);
                partial_syms = str2sym(partial_names);
                deps = obj.external_deps.(ext_names{i});
                for j = 1:length(deps)
                    % partial_deriv = functionalDerivative(obj.externals(i), deps(j));
                    partial_deriv = diff(exts(i), deps(j));
                    e = subs(e, partial_deriv, partial_syms(j));            
                end
                ext_d{i} = partial_syms;
            end
            ext_vars = str2sym(obj.getExternalName([], naming));
            e = subs(e, exts, ext_vars);

            [e, vars] = replaceDOFs(obj, e, naming);
            
            vars.ext = ext_vars;
            vars.ext_d = ext_d;

            vars.u = str2sym(obj.getInName([], naming));
            vars.p = str2sym(obj.getParamName([], naming));

            e = subs(e, struct2array(obj.inputs), vars.u);
            e = subs(e, struct2array(obj.params), vars.p);
        end

        % TODO: refactor to use one common jacobian function called with
        % different derivative vectors
        % Generalized mass matrix
        function M_ = getM(obj)
            if ~isempty(obj.M)
                M_ = obj.M;
                return
            end
            eom_ = obj.getEOM();
            M_= jacobian(eom_, diff(obj.q, obj.time, 2));
            obj.M = M_;
        end

        % Generalized coriolis and damping matrix
        function C_ = getC(obj)
            if ~isempty(obj.C)
                C_ = obj.C;
                return
            end
            eom_ = obj.getEOM();
            C_= jacobian(eom_, diff(obj.q, obj.time));
            obj.C = C_;
        end

        % Generalized stiffness matrix
        function K_ = getK(obj)
            if ~isempty(obj.K)
                K_ = obj.K;
                return
            end
            eom_ = obj.getEOM();
            K_= jacobian(eom_, obj.q);
            obj.K = K_;
        end

        % Input Jacobian
        function B_ = getB(obj)
            if ~isempty(obj.B)
                B_ = obj.B;
                return
            end
            eom_ = obj.getEOM();
            B_= jacobian(eom_, struct2array(obj.inputs));
            obj.B = B_;
        end

        % Output Jacobian
        function CD_ = getCD(obj)
            if ~isempty(obj.CD) || isempty(fieldnames(obj.outputs))
                CD_ = obj.CD;
                return
            end
            CD_= jacobian(struct2array(obj.outputs), [obj.q; diff(obj.q, obj.time); struct2array(obj.inputs)]);
            obj.CD = CD_;
        end

        % Output Jacobian wrt accelerations
        function F_ = getF(obj)
            if ~isempty(obj.F) || isempty(fieldnames(obj.outputs))
                F_ = obj.F;
                return
            end
            F_= jacobian(struct2array(obj.outputs), diff(obj.q, obj.time, 2));
            obj.F = F_;
        end

        function v = paramVec(obj, p)
            pv = struct2array(obj.params);
            v = zeros(size(pv));

            for i= 1:length(pv)
                v(i) = p.(char(pv(i)));
            end
        end
    end
    methods (Access = private)
        % calculate kinematics
        function prepareKinematics(obj)
            for i= 1:length(obj.children)
                obj.children(i).T0= obj.children(i).T;
                obj.children(i).prepareKinematics;
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
            if isfield(obj.params, name)
                if nargout>0
                    name_in_use= true;
                else
                    error("Parameter name '%s' already exists in the system.", name);
                end
            end
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

        function checkSetupCompleted(obj)
            if obj.setupCompleted
                error('Model setup has been completed. No further data can be added.')

            end
        end
    end
    methods (Static)
        function applyForce(F, b1)
            arguments
                F (3,1) sym
                b1 (1,1) Body
            end
            b1.F_ext = b1.F_ext + F;
        end

        function applyMoment(M, b1)
            arguments
                M (3,1) sym
                b1 (1,1) Body
            end
            b1.M_ext = b1.M_ext + M;
        end

        function forceBetween(F, b1, b2)
            arguments
                F (3,1) sym
                b1 (1,1) Body
                b2 (1,1) Body
            end
            MultiBodySystem.applyForce(F, b1)
            MultiBodySystem.applyForce(-F, b2)
        end

        function momentBetween(M, b1, b2)
            arguments
                M (3,1) sym
                b1 (1,1) Body
                b2 (1,1) Body
            end
            MultiBodySystem.applyMoment(M, b1)
            MultiBodySystem.applyMoment(-M, b2)
        end

        function applyElasticForce(Fe, b1)
            arguments
                Fe (:,1) sym
                b1 (1,1) ElasticBody
            end

            b1.Fe_ext = b1.Fe_ext + Fe;
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
