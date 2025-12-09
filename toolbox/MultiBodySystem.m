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
        outputs struct = struct()
        params struct = struct()
        children (1,:) Body = Body.empty

        q_ (:,1) sym = []
        qd_ (:,1) sym = []
        qdd_ (:,1) sym = []

        aux_state struct = struct()             % auxiliary state names
        aux_impl_ode (:,1) sym = []             % auxiliary implicit first order ode        
    end

    properties (Access = private)
    end

    methods
        % Constructor with optional name
        function obj = MultiBodySystem(name)
            if nargin > 0
                obj.Name = string(name);
            end
        end

        % Add a generalized coordinate by name
        function addGeneralizedCoordinate(obj, coordName)
            arguments
                obj
                coordName (1,:) char
            end
    
            % Check for duplicate name
            checkName(obj, coordName)
    
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
                coordName (1,:) char
            end
    
            % Check for duplicate name
            checkName(obj, coordName)
    
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
        function addParameter(obj, paramName)
            arguments
                obj
                paramName (1,:) char
            end
    
            % Check for duplicate name
            checkName(obj, paramName)
    
            % Define symbolic variable dynamically
            symVar = sym(paramName, 'real');
    
            % Store it
            obj.params.(paramName) = symVar;
        end

        % Add input by name
        function addInput(obj, inName)
            arguments
                obj
                inName (1,:) char
            end
    
            % Check for duplicate name
            checkName(obj, inName)
    
            % Define symbolic variable dynamically
            symVar = sym(inName, 'real');
    
            % Store it
            obj.inputs.(inName) = symVar;
        end

        % Add auxilliary state by name
        function addAuxState(obj, auxName)
            arguments
                obj
                auxName (1,:) char
            end
    
            % Check for duplicate name
            checkName(obj, auxName)
    
            % Define symbolic variable dynamically
            symVar = sym(auxName, 'real');
    
            % Store it
            obj.aux_state.(auxName) = symVar;
        end

        % Get the time derivative of auxilliary state
        function dx = getAuxStateD(obj, auxName)
            arguments
                obj
                auxName (1,:) char
            end

            dx = diff(obj.aux_state.(auxName), obj.time);
        end

        % Add auxilliary state by name
        function addAuxImplODE(obj, ode)
            obj.aux_impl_ode(end+1) = ode;
        end

        % Add a body with a unique name
        function addBody(obj, body)
            arguments
                obj
                body (1,1) Body
            end
    
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
                checkName(obj, body.Name)
            end

            obj.bodies.(body.Name) = body;
        end

        % Add a direct child
        function addChild(obj, body)
            arguments
                obj
                body (1,1) Body
            end

            addBody(obj, body);
            
            body.parent= obj;
            body.system= obj;
            obj.children(end+1)= body;
        end

        % Get number of bodies
        function n = numBodies(obj)
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
        
        function n = n_dof(obj)
            n = length(obj.q);
        end

        % calculate kinematics
        function prepareKinematics(obj)
            for i= 1:length(obj.children)
                obj.children(i).T0= obj.children(i).T;
                obj.children(i).prepareKinematics;
            end
        end

        function eom = getEOM(obj, replace_diff)
            if ~exist('replace_diff', 'var')
                replace_diff = true;
            end

            eom = sym(zeros(length(obj.q), 1));
            for i= 1:length(obj.children)
                eom = eom + obj.children(i).collectGenForces;
            end
            
            eom = simplify(eom);
            if replace_diff
                % TODO: do this in model finalization step
                if isempty(obj.qdd_)
                    obj.qdd_ = sym('qdd_', [length(obj.q), 1]);
                    obj.qd_ = sym('qd_', [length(obj.q), 1]);
                    obj.q_ = sym('q_', [length(obj.q), 1]);
                end

                eom = subs(eom, diff(obj.q, obj.time, 2), obj.qdd_);
                eom = subs(eom, diff(obj.q, obj.time), obj.qd_);
                eom = subs(eom, obj.q, obj.q_);
            end
        end

        function Fz = getConstraintForces(obj, replace_diff)
            if ~exist('replace_diff', 'var')
                replace_diff = true;
            end

            Fz = sym(zeros(length(obj.z), 1));
            for i= 1:length(obj.children)
                Fz = Fz + obj.children(i).collectConstrForces;
            end
            
            Fz = simplify(Fz);
            if replace_diff
                % TODO: do this in model finalization step
                if isempty(obj.qdd_)
                    obj.qdd_ = sym('qdd_', [length(obj.q), 1]);
                    obj.qd_ = sym('qd_', [length(obj.q), 1]);
                    obj.q_ = sym('q_', [length(obj.q), 1]);
                end

                Fz = subs(Fz, diff(obj.q, obj.time, 2), obj.qdd_);
                Fz = subs(Fz, diff(obj.q, obj.time), obj.qd_);
                Fz = subs(Fz, obj.q, obj.q_);
            end
        end

        function fun = eomFunO2(obj, filename)
            eom = getEOM(obj);
            if ~exist('filename', 'var')
                fun = matlabFunction(eom, 'Vars', {obj.q_, obj.qd_, obj.qdd_, struct2array(obj.inputs), struct2array(obj.params)});            
            else
                matlabFunction(eom, 'File', filename, 'Vars', {obj.q_, obj.qd_, obj.qdd_, struct2array(obj.inputs), struct2array(obj.params)});
            end
        end

        function fun = eomFunDescriptor(obj, filename)
            eom = getEOM(obj);
            [f_, M] = getF_(obj, eom);

            if ~exist('filename', 'var')
                fun = matlabFunction(M, f_, 'Vars', {obj.q_, obj.qd_, struct2array(obj.inputs), struct2array(obj.params)});            
            else
                matlabFunction(M, f_, 'File', filename, 'Vars', {obj.q_, obj.qd_, struct2array(obj.inputs), struct2array(obj.params)});
            end
        end

        function fun = eomDae(obj, filename)
            eom = getEOM(obj, false);
            

            for i = 1:length(obj.q)
                x1(i, 1) = symfun(str2sym(sprintf('x1_%d(time)', i)), obj.time);
                x2(i, 1) = symfun(str2sym(sprintf('x2_%d(time)', i)), obj.time);
            end
            x = [x1 ; x2];

            eom = subs(eom, diff(obj.q, obj.time, 2), diff(x2, obj.time));
            eom = subs(eom, diff(obj.q, obj.time), x2);
            eom = subs(eom, obj.q, x1);

            f_impl = [diff(x1, obj.time) == x2 ; eom; obj.aux_impl_ode];

            if ~exist('filename', 'var')
                fun = daeFunction(f_impl, x, in, struct2array(obj.params));
            else
                daeFunction(f_impl, x, struct2array(obj.inputs), struct2array(obj.params), 'File', filename);
            end
        end

        function M = getM(obj, eom)
            M= -jacobian(eom, obj.qdd_);
        end

        function [f_, M] = getF_(obj, eom, M)
            if ~exist('M', 'var')
                M = getM(obj, eom);
            end
            f_= simplify(eom + M*obj.qdd_);
        end

        function v = paramVec(obj, p)
            pv = struct2array(obj.params);
            v = zeros(size(pv));

            for i= 1:length(pv)
                v(i) = p.(char(pv(i)));
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
            if isfield(obj.aux_state, name)
                if nargout>0
                    name_in_use= true;
                else
                    error("Auxilliary state name '%s' already exists in the system.", name);
                end
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
end
