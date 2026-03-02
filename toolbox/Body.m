classdef Body  < handle & matlab.mixin.Heterogeneous
    properties
        Name (1,:) char = ''
        Description (1,1) string = "No description"

        T                             % Homogeneous transformation matrix relative to parent. Multiply body coordinates with this matrix to get global coordinates
        children (1,:) Body = Body.empty         
        parent = []
        system (1,:) MultiBodySystem = MultiBodySystem.empty

        T0 (:,:)                      % Homogeneous transformation matrix from inertial reference frame
        v0 (3,1)
        v0_z (3,1)                    % Velocity including movement in constraint direction
        a0 (3,1) 
        omega0 (3,1)
        omega0_z (3,1)                % Rotational velocity including movement in constraint direction
        alpha0 (3,1) 
        v_p (3,:)
        omega_p (3,:)
        vz_p (3,:)                    % Partial velocities in constrain direction
        omegaz_p (3,:)                % Partial rotational velocities in constrain direction

        F_ext (:, 1) 
        M_ext (:, 1)
        
        F (:, 1) 
        M (:, 1)

        Fgen (:, 1) 
        Fconstr (:, 1)                % Constraint forces

        forcesPrepared = false        % replace this by a general locking/finalization of the entire system
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
            % Initialize with numeric matrices that will be converted to symbolic by the system
            obj.T = eye(4);
            obj.T0 = eye(4);
            obj.F_ext = [0 0 0].';
            obj.M_ext = [0 0 0].';
            obj.F = [0 0 0].';
            obj.M = [0 0 0].';
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
            T_translate = obj.system.sym(eye(4));
            T_translate(1:3,4) = vec;
            obj.T = obj.T * T_translate;
        end

        % Rotate about local X, Y, or Z axis by angle (in radians)
        function rotateLocalAxis(obj, axis, angle)
            axis = upper(axis);
            switch axis
                case 'X'
                    R = obj.rotationMatrix([1 0 0], angle);
                case 'Y'
                    R = obj.rotationMatrix([0 1 0], angle);
                case 'Z'
                    R = obj.rotationMatrix([0 0 1], angle);
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
            R = obj.rotationMatrix(axisNorm, angle);
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
            obj.system.checkKinematicsFinished()

            r_abs = obj.T0(1:3, 1:3) * r_rel;
            v_abs = obj.T0(1:3, 1:3) * v_rel;
            a_abs = obj.T0(1:3, 1:3) * a_rel;
        
            abs_accel = obj.system.simplify(obj.a0 + crossmat(obj.alpha0)*r_abs + crossmat(obj.omega0)*(crossmat(obj.omega0)*r_abs) + 2*crossmat(obj.omega0)*v_abs + a_abs);
        end

        % calculate kinematics - dispatcher method
        function prepareKinematicsBase(obj, kinematics_from_global)
            if kinematics_from_global
                obj.prepareKinematicsFromGlobal();
            else
                obj.prepareKinematicsFromLocal();
            end

            for i= 1:length(obj.children)
                obj.children(i).T0= obj.T0 * obj.children(i).T;
                obj.children(i).prepareKinematics(kinematics_from_global);
            end

            obj.T0= obj.system.removeDOC(obj.T0);

            % eps may only be removed after all children have been processed, otherwise recursive propagation is not correct.
            if ~kinematics_from_global
                obj.v_p = obj.system.simplify(obj.system.removeEps(obj.v_p));
                obj.omega_p = obj.system.simplify(obj.system.removeEps(obj.omega_p));
                obj.vz_p = obj.system.simplify(obj.system.removeEps(obj.vz_p));
                obj.omegaz_p = obj.system.simplify(obj.system.removeEps(obj.omegaz_p));
            end
        end

        % calculate kinematics from global frame
        function prepareKinematicsFromGlobal(obj)
            obj.v0_z= obj.system.simplify(obj.system.removeEps(diff(obj.T0(1:3, 4), obj.system.time), true));
            obj.v0 = obj.system.removeDOC(obj.v0_z);

            obj.a0= diff(obj.v0, obj.system.time);

            w_skew= obj.system.simplify(diff(obj.T0(1:3, 1:3), obj.system.time)*obj.T0(1:3, 1:3).');
            obj.omega0_z= [(w_skew(3, 2)-w_skew(2, 3))/obj.system.sym(2); (w_skew(1, 3)-w_skew(3, 1))/obj.system.sym(2); (w_skew(2, 1)-w_skew(1, 2))/obj.system.sym(2)];
            obj.omega0_z= obj.system.simplify(obj.system.removeEps(obj.omega0_z, true));
            obj.omega0 = obj.system.removeDOC(obj.omega0_z);

            obj.alpha0= diff(obj.omega0, obj.system.time);

            % Partial velocities
            % TODO: add explanation, why eps is removed here
            obj.v_p= obj.system.simplify(obj.system.keepEps(jacobian(obj.v0, diff(obj.system.q, obj.system.time))));
            obj.omega_p= obj.system.simplify(obj.system.keepEps(jacobian(obj.omega0, diff(obj.system.q, obj.system.time))));

            obj.vz_p= obj.system.keepEps(jacobian(obj.v0_z, diff(obj.system.z, obj.system.time)));
            obj.vz_p= obj.system.simplify(obj.system.removeDOC(obj.vz_p));

            obj.omegaz_p= obj.system.keepEps(jacobian(obj.omega0_z, diff(obj.system.z, obj.system.time)));
            obj.omegaz_p= obj.system.simplify(obj.system.removeDOC(obj.omegaz_p));
        end

        % calculate kinematics from local frame
        function prepareKinematicsFromLocal(obj)
            % Calculate local kinematics first
            % Velocity and acceleration from local transformation
            v_local_z = diff(obj.T(1:3, 4), obj.system.time);
            v_local = obj.system.removeDOC(obj.system.removeEps(v_local_z, true));
            a_local = diff(v_local, obj.system.time);
            
            % Angular velocity from local rotation matrix
            w_skew_local = obj.system.simplify(diff(obj.T(1:3, 1:3), obj.system.time) * obj.T(1:3, 1:3).');
            omega_local_z = [(w_skew_local(3, 2)-w_skew_local(2, 3))/obj.system.sym(2); ...
                          (w_skew_local(1, 3)-w_skew_local(3, 1))/obj.system.sym(2); ...
                          (w_skew_local(2, 1)-w_skew_local(1, 2))/obj.system.sym(2)];
            omega_local_z = obj.system.simplify(obj.system.removeEps(omega_local_z, true));
            omega_local = obj.system.removeDOC(omega_local_z);
            alpha_local = diff(omega_local, obj.system.time);
            
            % Partial velocities
            v_local_p = jacobian(v_local, diff(obj.system.q, obj.system.time));
            omega_local_p = jacobian(omega_local, diff(obj.system.q, obj.system.time));

            v_local_z_p = obj.system.removeDOC(jacobian(v_local_z, diff(obj.system.z, obj.system.time)));
            omega_local_z_p= obj.system.removeDOC(jacobian(omega_local_z, diff(obj.system.z, obj.system.time)));

            % omega0_z and v0_z are not needed, so leave them out
            if isa(obj.parent, 'MultiBodySystem')
                % If parent is MultiBodySystem = inertial frame, local and global kinematics are the same
                obj.v0 = v_local;
                obj.a0 = a_local;
                obj.omega0 = omega_local;
                obj.alpha0 = alpha_local;

                obj.v_p = v_local_p;
                obj.omega_p = omega_local_p;
                obj.vz_p = v_local_z_p;
                obj.omegaz_p = omega_local_z_p;
            else
                rotationToGlobal = obj.system.removeDOC(obj.parent.T0(1:3, 1:3));
                r_rel = rotationToGlobal * obj.system.removeDOC(obj.T(1:3, 4));
                v_rel = rotationToGlobal * v_local;
                omega_rel = rotationToGlobal * omega_local;
                a_rel = rotationToGlobal * a_local;
                alpha_rel = rotationToGlobal * alpha_local;

                obj.v0 = obj.parent.v0 + crossmat(obj.parent.omega0) * r_rel + v_rel;
                obj.a0 = obj.parent.a0 + crossmat(obj.parent.alpha0) * r_rel + crossmat(obj.parent.omega0) * (crossmat(obj.parent.omega0) * r_rel) + 2 * crossmat(obj.parent.omega0) * v_rel + a_rel;
                obj.omega0 = obj.parent.omega0 + omega_rel;
                obj.alpha0 = obj.parent.alpha0 + crossmat(obj.parent.omega0) * omega_rel + alpha_rel;

                obj.v_p = obj.system.sym(zeros(3, length(obj.system.q)));
                obj.omega_p = obj.system.sym(zeros(3, length(obj.system.q)));
                for i= 1:length(obj.system.q)
                    v_rel_p = rotationToGlobal * v_local_p(:, i);
                    obj.v_p(:, i) = obj.parent.v_p(:, i) + crossmat(obj.parent.omega_p(:, i)) * r_rel + v_rel_p;

                    omega_rel_p = rotationToGlobal * omega_local_p(:, i);
                    obj.omega_p(:, i) = obj.parent.omega_p(:, i) + omega_rel_p;
                end
                obj.vz_p = obj.system.sym(zeros(3, length(obj.system.z)));
                obj.omegaz_p = obj.system.sym(zeros(3, length(obj.system.z)));
                for i= 1:length(obj.system.z)
                    v_rel_z_p = rotationToGlobal * v_local_z_p(:, i);
                    obj.vz_p(:, i) = obj.parent.vz_p(:, i) + crossmat(obj.parent.omegaz_p(:, i)) * r_rel + v_rel_z_p;

                    omega_rel_z_p = rotationToGlobal * omega_local_z_p(:, i);
                    obj.omegaz_p(:, i) = obj.parent.omegaz_p(:, i) + omega_rel_z_p;
                end
            end
        end

        function prepareForces(~)
            % virtual method implemented by subclasses
        end

        function calcGenForce(obj)
            obj.Fgen = - obj.v_p.' * (obj.F + obj.F_ext);
            obj.Fgen = obj.Fgen - obj.omega_p.' * (obj.M + obj.M_ext);
            obj.Fgen = obj.system.simplify(obj.Fgen);
        end

        function calcConstrForce(obj)
            obj.Fconstr = obj.vz_p.' * (obj.F + obj.F_ext);
            obj.Fconstr = obj.Fconstr + obj.omegaz_p.' * (obj.M + obj.M_ext);
            obj.Fconstr = obj.system.simplify(obj.Fconstr);
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
                F (3,1) 
            end
            obj.system.checkKineticsNotFinished();

            obj.F_ext = obj.F_ext + F;
        end

        % Apply a moment given in global coordinates to this body
        function applyMoment(obj, M)
            arguments
                obj
                M (3,1) 
            end
            obj.system.checkKineticsNotFinished();

            obj.M_ext = obj.M_ext + M;
        end

        % Apply a force given in global coordinates to this body and
        % negative to another body
        function forceBetween(obj, F, b2)
            arguments
                obj
                F (3,1) 
                b2 (1,1) Body
            end
            obj.system.checkKineticsNotFinished();

            obj.applyForce(F)
            b2.applyForce(-F)
        end

        % Apply a moment given in global coordinates to this body and
        % negative to another body
        function momentBetween(obj, M, b2)
            arguments
                obj
                M (3,1) 
                b2 (1,1) Body
            end
            obj.system.checkKineticsNotFinished();

            obj.applyMoment(M)
            b2.applyMoment(-M)
        end
        
        % Apply a force given in body coordinates at a point on the body
        % given in body coordinates
        function applyForceInLocal(obj, r, F)
            arguments
                obj
                r (3,1)          % position relative to center of mass or reference system in body local coordinates
                F (3,1)          % force in body local coordinates
            end
            % make sure T0 is already available
            obj.system.checkKinematicsFinished()
            obj.system.checkKineticsNotFinished()

            Fin0 = obj.T0(1:3, 1:3) * F;
            r0 = obj.T0(1:3, 1:3) * r;
            obj.applyForceIn0(r0, Fin0)
        end
    
        % Apply a force given in global coordinates at a point on the body
        % given in global coordinates
        function applyForceIn0(obj, r0, F0)
            arguments
                obj
                r0 (3,1)          % position relative to center of mass or reference system in global coordinates
                F0 (3,1)          % force in global coordinates
            end
            obj.system.checkKineticsNotFinished()

            obj.applyForce(F0)
            obj.applyMoment(crossmat(r0)*F0)
        end

        % Rodrigues' formula for rotation matrix and wrap into 4x4
        % TODO: move to MultiBodySystem
        function T = rotationMatrix(obj, axis, angle)
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

            T = obj.system.sym(eye(4));
            T(1:3,1:3) = R;
        end
    end

    methods (Static)
        function mustBeNumericOrSym(val)
            if ~(isnumeric(val) || isa(val,'msym') || isa(val,'sym'))
                error('Value must be numeric or symbolic.');
            end
        end
    end
end
