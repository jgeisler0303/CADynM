classdef Body  < handle & matlab.mixin.Heterogeneous
    properties
        Name (1,:) char = ''
        Description (1,1) string = "No description"

        T (:,:) msym = []                        % Homogeneous transformation matrix relative to parent. Multiply body coordinates with this matrix to get global coordinates
        children (1,:) Body = Body.empty         
        parent = []
        system (1,:) MultiBodySystem = MultiBodySystem.empty

        T0 (:,:) msym = []                   % Homogeneous transformation matrix from inertial reference frame
        v0 (3,1) msym = [0 0 0]'
        v0_z (3,1) msym = [0 0 0]'               % Velocity including movement in constraint direction
        a0 (3,1) msym = [0 0 0]'
        omega0 (3,1) msym = [0 0 0]'
        omega0_z (3,1) msym = [0 0 0]'           % Rotational velocity including movement in constraint direction
        alpha0 (3,1) msym = [0 0 0]'
        v_p (3,:) msym = msym.empty(3, 0)
        omega_p (3,:) msym = msym.empty(3, 0)
        vz_p (3,:) msym = msym.empty(3, 0)       % Partial velocities in constrain direction
        omegaz_p (3,:) msym = msym.empty(3, 0)   % Partial rotational velocities in constrain direction

        F_ext (:, 1) msym = []
        M_ext (:, 1) msym = []
        
        F (:, 1) msym = []
        M (:, 1) msym = []

        Fgen (:, 1) msym = []
        Fconstr (:, 1) msym = []                % Constraint forces

        forcesPrepared = false                  % replace this by a general locking/finalization of the entire system
    end

    methods
        % Constructor
        function obj = Body(name, description)
            if nargin > 0
                obj.Name = name;
            end
            if nargin > 1
                obj.Description = description;
            end
            obj.T = eye(4);
            obj.T0 = eye(4);
            obj.F_ext = [0 0 0]';
            obj.M_ext = [0 0 0]';
            obj.F = [0 0 0]';
            obj.M = [0 0 0]';
        end

        % Add a direct child
        function addChild(obj, body)
            arguments
                obj
                body (1,1) Body
            end

            obj.system.addBody(body);
            
            body.parent= obj;
            body.system= obj.system;
            obj.children(end+1)= body;
        end

        % Apply translation (vector can be row or column)
        function translate(obj, vec)
            vec = vec(:);  % ensure column vector
            if numel(vec) ~= 3
                error("Translation vector must have 3 elements.");
            end
            T_translate = msym(eye(4));
            T_translate(1:3,4) = vec;
            obj.T = obj.T * T_translate;
        end

        % Rotate about local X, Y, or Z axis by angle (in radians)
        function rotateLocalAxis(obj, axis, angle)
            axis = upper(axis);
            switch axis
                case 'X'
                    R = Body.rotationMatrix([1 0 0], angle);
                case 'Y'
                    R = Body.rotationMatrix([0 1 0], angle);
                case 'Z'
                    R = Body.rotationMatrix([0 0 1], angle);
                otherwise
                    error("Axis must be 'X', 'Y', or 'Z'.");
            end
            obj.T = obj.T * R;
        end

        % Rotate about an arbitrary axis in the current body frame
        function rotateAboutAxis(obj, axisVec, angle)
            axisVec = axisVec(:);
            if numel(axisVec) ~= 3 || norm(axisVec) == 0
                error("Rotation axis must be a non-zero 3D vector.");
            end
            axisNorm = axisVec / norm(axisVec);
            R = Body.rotationMatrix(axisNorm, angle);
            obj.T = obj.T * R;
        end

        % calculate acceleration of local point in global coordinates
        % TODO: extend to any reference frame
        function abs_accel = getA0(obj, r_rel, v_rel, a_rel)
            arguments
                obj 
                r_rel (3, 1) {Body.mustBeNumericOrSym} 
                v_rel (3, 1) {Body.mustBeNumericOrSym} 
                a_rel (3, 1) {Body.mustBeNumericOrSym} 
            end
            obj.system.checkSetupCompleted()

            r_abs = obj.T0(1:3, 1:3) * r_rel;
            v_abs = obj.T0(1:3, 1:3) * v_rel;
            a_abs = obj.T0(1:3, 1:3) * a_rel;
        
            abs_accel = simplify(obj.a0 + crossmat(obj.alpha0)*r_abs + crossmat(obj.omega0)*(crossmat(obj.omega0)*r_abs) + 2*crossmat(obj.omega0)*v_abs + a_abs);
        end

        % calculate kinematics
        function prepareKinematicsBase(obj)
            obj.v0= diff(obj.T0(1:3, 4), obj.system.time);
            obj.v0_z = obj.v0;
            % store with and remove movement inconstraint directions
            obj.v0= subs(obj.v0, struct2array(obj.system.doc), 0);

            obj.a0= diff(obj.v0, obj.system.time);

            w_skew= simplify(diff(obj.T0(1:3, 1:3), obj.system.time)*obj.T0(1:3, 1:3).');
            obj.omega0= [(w_skew(3, 2)-w_skew(2, 3))/msym(2); (w_skew(1, 3)-w_skew(3, 1))/msym(2); (w_skew(2, 1)-w_skew(1, 2))/msym(2)];
            obj.omega0= simplify(obj.system.removeEps(obj.omega0, true));
            % store with and remove movement inconstraint directions
            obj.omega0_z= obj.omega0;
            obj.omega0 = subs(obj.omega0, struct2array(obj.system.doc), 0);

            obj.alpha0= diff(obj.omega0, obj.system.time);

            % Partial velocities
            % TODO: add explanation, why eps is removed here
            obj.v_p= simplify(subs(jacobian(obj.v0, diff(obj.system.q, obj.system.time)), [obj.system.sym_eps, obj.system.sym_eps_rot], [1, 1]));
            obj.omega_p= simplify(subs(jacobian(obj.omega0, diff(obj.system.q, obj.system.time)), [obj.system.sym_eps, obj.system.sym_eps_rot], [1, 1]));

            obj.vz_p= subs(jacobian(obj.v0_z, diff(obj.system.z, obj.system.time)), [obj.system.sym_eps, obj.system.sym_eps_rot], [1, 1]);
            obj.vz_p= simplify(subs(obj.vz_p, struct2array(obj.system.doc), 0));

            obj.omegaz_p= subs(jacobian(obj.omega0_z, diff(obj.system.z, obj.system.time)), [obj.system.sym_eps, obj.system.sym_eps_rot], [1, 1]);
            obj.omegaz_p= simplify(subs(obj.omegaz_p, struct2array(obj.system.doc), 0));

            for i= 1:length(obj.children)
                obj.children(i).T0= obj.T0 * obj.children(i).T;
                obj.children(i).prepareKinematics;
            end

            obj.T0= subs(obj.T0, struct2array(obj.system.doc), 0);
        end

        function prepareForces(obj)
        end

        function calcGenForce(obj)
            obj.Fgen = - obj.v_p.' * (obj.F + obj.F_ext);
            obj.Fgen = obj.Fgen - obj.omega_p.' * (obj.M + obj.M_ext);
            obj.Fgen = simplify(obj.Fgen);
        end

        function calcConstrForce(obj)
            obj.Fconstr = obj.vz_p.' * (obj.F + obj.F_ext);
            obj.Fconstr = obj.Fconstr + obj.omegaz_p.' * (obj.M + obj.M_ext);
            obj.Fconstr = simplify(obj.Fconstr);
        end

        function Fgen = collectGenForces(obj)
            if ~obj.forcesPrepared
                obj.prepareForces();
                obj.forcesPrepared = true;
            end
            obj.calcGenForce();
            
            obj.Fgen = obj.system.removeEps(obj.Fgen);
            
            Fgen = obj.Fgen;
            for i= 1:length(obj.children)
                Fgen = Fgen + obj.children(i).collectGenForces;
            end
        end

        function Fconstr = collectConstrForces(obj)
            if ~obj.forcesPrepared
                obj.prepareForces();
                obj.forcesPrepared = true;
            end
            obj.calcConstrForce();

            Fconstr = obj.Fconstr;
            for i= 1:length(obj.children)
                Fconstr = Fconstr + obj.children(i).collectConstrForces;
            end
        end

        % Apply a force given in global coordinates to the center of gravitiy of this body
        function applyForce(obj, F)
            arguments
                obj
                F (3,1) msym
            end
            obj.F_ext = obj.F_ext + F;
        end

        % Apply a moment given in global coordinates to this body
        function applyMoment(obj, M)
            arguments
                obj
                M (3,1) msym
            end
            obj.M_ext = obj.M_ext + M;
        end

        % Apply a force given in global coordinates to this body and
        % negative to another body
        function forceBetween(obj, F, b2)
            arguments
                obj
                F (3,1) msym
                b2 (1,1) Body
            end
            obj.applyForce(F)
            b2.applyForce(-F)
        end

        % Apply a moment given in global coordinates to this body and
        % negative to another body
        function momentBetween(obj, M, b2)
            arguments
                obj
                M (3,1) msym
                b2 (1,1) Body
            end
            obj.applyMoment(M)
            b2.applyMoment(-M)
        end
        
        % Apply a force given in body coordinates at a point on the body
        % given in body coordinates
        function applyForceInLocal(obj, r, F)
            arguments
                obj
                r (3,1) msym         % position relative to center of mass or reference system in body local coordinates
                F (3,1) msym         % force in body local coordinates
            end
            % make sure T0 is already available
            obj.system.checkSetupCompleted()

            Fin0 = obj.T0(1:3, 1:3) * F;
            r0 = obj.T0(1:3, 1:3) * r;
            obj.applyForceIn0(r0, Fin0)
        end
    
        % Apply a force given in global coordinates at a point on the body
        % given in global coordinates
        function applyForceIn0(obj, r0, F0)
            arguments
                obj
                r0 (3,1) msym         % position relative to center of mass or reference system in global coordinates
                F0 (3,1) msym         % force in global coordinates
            end
            obj.applyForce(F0)
            obj.applyMoment(crossmat(r0)*F0)
        end        
    end

    methods (Static)
        % Rodrigues' formula for rotation matrix and wrap into 4x4
        function T = rotationMatrix(axis, angle)
            if ischar(axis)
                switch axis
                    case 'x'
                        axis = [1 0 0];
                    case 'y'
                        axis = [0 1 0];
                    case 'z'
                        axis = [0 0 1];
                    otherwise
                        error('Unknown rotation axis name "%s".', axis)
                end
            end
            x = axis(1); y = axis(2); z = axis(3);
            c = cos(angle);
            s = sin(angle);
            C = 1 - c;

            R = [x*x*C + c,     x*y*C - z*s, x*z*C + y*s;
                 y*x*C + z*s, y*y*C + c,     y*z*C - x*s;
                 z*x*C - y*s, z*y*C + x*s, z*z*C + c];

            T = msym(eye(4));
            T(1:3,1:3) = R;
        end
        
        function mustBeNumericOrSym(val)
            if ~(isnumeric(val) || isa(val,'msym'))
                error('Value must be numeric or symbolic.');
            end
        end
    end
end
