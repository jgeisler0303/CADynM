classdef RigidBody  < Body
    properties
        m (1,1) sym = 0                         % Mass
        I (3,3) sym = zeros(3,3)                % Inertia matrix
    end

    methods
        % Constructor
        function obj = RigidBody(name, description, mass, inertia)
            if nargin > 0
                obj.Name = string(name);
            end
            if nargin > 1
                if ~isempty(description)
                    obj.Description = string(description);
                end
            end
            if nargin > 2
                obj.m = mass;
            end
            if nargin > 3
                obj.I = inertia;
            end
        end

        function prepareKinematics(obj)
            prepareKinematicsBase(obj)
        end

        function prepareForces(obj)
            obj.F = obj.m * (obj.system.gravity - obj.a0);

            rotationToGlobal = obj.T0(1:3,1:3);
            rotationToLocal = rotationToGlobal.';

            Phi = obj.I;
            PhiG_global = rotationToGlobal * Phi * rotationToLocal;

            obj.M = -PhiG_global * obj.alpha0 - crossmat(obj.omega0) * (PhiG_global * obj.omega0);
        end
    end
end
