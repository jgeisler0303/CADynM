classdef Body  < handle & matlab.mixin.Heterogeneous
    properties
        Name (1,:) char = ''
        Description (1,1) string = "No description"

        T (4,4) sym = eye(4)                    % Homogeneous transformation matrix
        children (1,:) Body = Body.empty         
        parent = []
        system (1,:) MultiBodySystem = MultiBodySystem.empty

        T0 (4,4) sym = eye(4)                    % Homogeneous transformation matrix from inertial reference frame
        v0 (3,1) sym = 0
        a0 (3,1) sym = 0
        omega0 (3,1) sym = 0
        alpha0 (3,1) sym = 0 
        v_p (3,:) sym = []
        omega_p (3,:) sym = []

        F_ext (3, 1) sym = 0
        M_ext (3, 1) sym = 0
        
        F (3, 1) sym = 0
        M (3, 1) sym = 0

        Fgen (:, 1) sym = 0
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
            T_translate = sym(eye(4));
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

        % calculate kinematics
        function prepareKinematicsBase(obj)
            obj.v0= diff(obj.T0(1:3, 4), obj.system.time);
            obj.a0= diff(obj.v0, obj.system.time);

            w_skew= simplify(diff(obj.T0(1:3, 1:3), obj.system.time)*obj.T0(1:3, 1:3).');
            obj.omega0= [[0 0 1]*w_skew*[0 1 0]'; [1 0 0]*w_skew*[0 0 1]'; [0 1 0]*w_skew*[1 0 0]'];
            obj.alpha0= diff(obj.omega0, obj.system.time);

            obj.v_p= jacobian(obj.v0, diff(obj.system.q, obj.system.time));
            obj.omega_p= jacobian(obj.omega0, diff(obj.system.q, obj.system.time));

            for i= 1:length(obj.children)
                obj.children(i).T0= obj.T0 * obj.children(i).T;
                obj.children(i).prepareKinematics;
            end
        end

        function prepareForces(obj)
        end

        function calcGenForce(obj)
            obj.Fgen = obj.v_p.' * (obj.F + obj.F_ext);
            obj.Fgen = obj.Fgen + obj.omega_p.' * (obj.M + obj.M_ext);            
        end

        function Fgen = collectGenForces(obj)
            obj.prepareForces();
            obj.calcGenForce();


            Fgen = expand(obj.Fgen);
            Fgen = mapSymType(Fgen, 'power', @(Z) subs(Z, sym('eps'), 0));
            Fgen = subs(Fgen, sym('eps'), 1);
            Fgen = simplify(Fgen);

            for i= 1:length(obj.children)
                Fgen = Fgen + obj.children(i).collectGenForces;
            end
        end
    end

    methods (Static, Access = private)
        % Rodrigues' formula for rotation matrix and wrap into 4x4
        function T = rotationMatrix(axis, angle)
            x = axis(1); y = axis(2); z = axis(3);
            c = cos(angle);
            s = sin(angle);
            C = 1 - c;

            R = [x*x*C + c,     x*y*C - z*s, x*z*C + y*s;
                 y*x*C + z*s, y*y*C + c,     y*z*C - x*s;
                 z*x*C - y*s, z*y*C + x*s, z*z*C + c];

            T = sym(eye(4));
            T(1:3,1:3) = R;
        end
    end
end
