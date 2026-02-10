classdef ElasticBody  < Body
    properties
        ElasticDOF (:,1) msym = []
        sid SID
        Fe (:, 1) msym = []
        Fe_ext (:, 1) msym = []
        e_p (:,:) msym = []
        ConstrLoads (1, :) cell = {}
        ChildFrame (1,:) double = []    % number of frame to which a child is attached
    end

    methods
        % Constructor
        function obj = ElasticBody(eDOF, sid, name, description)
            if isstruct(sid)
                obj.sid = SID(sid);
            else
                obj.sid = sid;
            end

            if ~isempty(eDOF) && length(eDOF)~=obj.sid.nelastq
                error('Number of elastic DOF (%d) and DOF in SID (%d) must match.', length(eDOF), sid.nelastq)
            end
            obj.ElasticDOF = eDOF;

            if nargin > 2
                obj.Name = name;
            end
            if nargin > 3
                obj.Description = description;
            end
        end

        % attach a body to one of this bodies frames
        function addChild(obj, body, iframe, with_loads)
            arguments
                obj
                body (1,1) Body
                iframe = []
                with_loads (1,1) logical = true
            end
            if isempty(iframe) || (ischar(iframe) && strcmpi(iframe, 'last'))
                iframe = length(obj.sid.frame);
            end
            if with_loads
                % TODO: also consider load moments
                loadConstrNames = strcat('loadConstr_', body.Name, '_on_' , obj.Name, '_', {'x' 'y' 'z'});
                obj.system.addConstraintCoordinate(loadConstrNames);
                obj.ConstrLoads{end+1} = loadConstrNames;
                obj.ChildFrame(end+1) = iframe;
                z_load = [obj.system.doc.(loadConstrNames{1}); obj.system.doc.(loadConstrNames{2}); obj.system.doc.(loadConstrNames{3})];
            else
                obj.ConstrLoads{end+1} = {};
                obj.ChildFrame(end+1) = -1;
                z_load = [0 0 0]';
            end

            T_elast= msym(eye(4));
            T_elast(1:3, 4) = obj.sid.frame(iframe).origin.evalTaylor(obj.ElasticDOF, obj.system.sym_eps) + z_load;
            T_elast(1:3, 1:3) = Trot_elast(obj, iframe);
            body.T = T_elast * body.T; % pre-multiply because all body transformations come after the movement caused by being attached to the elastic body
        
            obj.system.addBody(body);
            body.parent= obj;
            body.system= obj.system;
            obj.children(end+1)= body;
        end

        % Apply an elastic force to the generalized coordinates of this
        % body
        % the force vector has to have length equal to the number of
        % elastic coordinates.
        % TODO: add function to apply force to node to be multiplied by phi
        function applyElasticForce(obj, Fe)
            arguments
                obj
                Fe (:,1) msym
            end

            obj.Fe_ext = obj.Fe_ext + Fe;
        end

        % Rotation matrix for elastic deformation
        function R = Trot_elast(obj, iframe)
            arguments
                obj
                iframe = []
            end
            if isempty(iframe) || (ischar(iframe) && strcmpi(iframe, 'last'))
                iframe = length(obj.sid.frame);
            end
            R = obj.sid.frame(iframe).ap.evalTaylor(obj.ElasticDOF, obj.system.sym_eps_rot*obj.system.sym_eps);
        end
        
        % calculate acceleration of elastic point in global coordinates
        % TODO: extend to any reference frame
        function abs_accel = getA0Elast(obj, iframe)
            arguments
                obj 
                iframe = []
            end
            if isempty(iframe) || (ischar(iframe) && strcmpi(iframe, 'last'))
                iframe = length(obj.sid.frame);
            end
            
            obj.system.checkSetupCompleted()

            r_rel = obj.sid.frame(iframe).origin.evalTaylor(obj.ElasticDOF, obj.system.sym_eps);
            v_rel = diff(r_rel, obj.system.time);
            a_rel = diff(v_rel, obj.system.time);
        
            abs_accel = obj.getA0(r_rel, v_rel, a_rel);
        end
        
        function prepareKinematics(obj)
            obj.e_p= jacobian(obj.ElasticDOF, obj.system.q);
            prepareKinematicsBase(obj)
        end

        function prepareForces(obj)
            eDOF_d= diff(obj.ElasticDOF, obj.system.time);
            eDOF_dd= diff(eDOF_d, obj.system.time);
            
            rotationToGlobal = obj.T0(1:3,1:3);
            rotationToLocal = rotationToGlobal.';
            
            obj.F = obj.sid.mass * (obj.system.gravity - obj.a0);
            
            md_global = rotationToGlobal * obj.sid.md.evalTaylor(obj.ElasticDOF, obj.system.sym_eps);
            obj.F = obj.F - crossmat(obj.alpha0) * md_global;
            obj.F = obj.F - crossmat(obj.omega0) * (crossmat(obj.omega0) * md_global);

            Ct = evalTaylor(obj.sid.Ct, obj.ElasticDOF, obj.system.sym_eps);
            Ct_d = zeros(3, 1);
            Ct_dd = zeros(3, 1);
            for ief = 1:length(obj.ElasticDOF)
                Ct_d = Ct_d + Ct(ief, :).' * eDOF_d(ief)*obj.system.sym_eps;
                Ct_dd = Ct_dd + Ct(ief, :).' * eDOF_dd(ief)*obj.system.sym_eps;
            end
            obj.F = obj.F - 2 * crossmat(obj.omega0) * rotationToGlobal * Ct_d;
            obj.F = obj.F - rotationToGlobal * Ct_dd;

            Phi = evalTaylor(obj.sid.I, obj.ElasticDOF, obj.system.sym_eps);

            PhiG_global = rotationToGlobal * Phi * rotationToLocal;

            obj.M = -PhiG_global * obj.alpha0 - crossmat(obj.omega0) * (PhiG_global * obj.omega0);

            omega_local = rotationToLocal * obj.omega0;

            obj.M = obj.M + crossmat(md_global) * (obj.system.gravity - obj.a0);

            Cr = evalTaylor(obj.sid.Cr, obj.ElasticDOF, obj.system.sym_eps);
            Cr_dd = zeros(3, 1);
            for ief = 1:length(obj.ElasticDOF)
                Cr_dd = Cr_dd + Cr(ief, :).' * eDOF_dd(ief)*obj.system.sym_eps;
            end
            obj.M = obj.M - rotationToGlobal * Cr_dd;

            Gr = evalTaylor(obj.sid.Gr, obj.ElasticDOF, obj.system.sym_eps);
            Gr_ = zeros(3, 3);
            for ief = 1:length(obj.ElasticDOF) 
                Gr_ = Gr_ + Gr{ief}(:, :) * eDOF_d(ief)*obj.system.sym_eps;
            end
            obj.M = obj.M - rotationToGlobal * Gr_ * omega_local;

            grav_local = rotationToLocal * (obj.a0 - obj.system.gravity);
            alpha_local = rotationToLocal * obj.alpha0;

            obj.Fe = zeros(length(obj.ElasticDOF), 1);

            for ief = 1:length(obj.ElasticDOF)
                obj.Fe(ief) = - Ct(ief, :) * grav_local;
                obj.Fe(ief) = obj.Fe(ief) - Cr(ief, :) * alpha_local;
                for jef = 1:length(obj.ElasticDOF)
                    obj.Fe(ief) = obj.Fe(ief) - obj.sid.Me.M0(ief, jef) * eDOF_dd(jef)*obj.system.sym_eps;
                end
            end

            w = [omega_local(1)*omega_local(1), omega_local(2)*omega_local(2), omega_local(3)*omega_local(3), omega_local(1)*omega_local(2), omega_local(2)*omega_local(3), omega_local(1)*omega_local(3)];
            obj.Fe = obj.Fe - evalTaylor(obj.sid.Oe, obj.ElasticDOF, obj.system.sym_eps) * w.';

            Ge = evalTaylor(obj.sid.Ge, obj.ElasticDOF, obj.system.sym_eps);
            for ief = 1:length(obj.ElasticDOF)
                Ge_ = zeros(1, 3);
                for jef = 1:length(obj.ElasticDOF)
                    % TODO = check for correct order of jef and ief
                    Ge_ = Ge_ + Ge{jef}(ief, :) * eDOF_d(jef)*obj.system.sym_eps;
                end
                obj.Fe(ief) = obj.Fe(ief) - Ge_ * omega_local;
            end

            K = evalTaylor(obj.sid.Ke, obj.ElasticDOF, obj.system.sym_eps);
            D = evalTaylor(obj.sid.De, obj.ElasticDOF, obj.system.sym_eps);
            for ief = 1:length(obj.ElasticDOF)
                for jef = 1 : length(obj.ElasticDOF)
                    obj.Fe(ief) = obj.Fe(ief) - K(ief, jef) * obj.ElasticDOF(jef)*obj.system.sym_eps;  % eps?
                    obj.Fe(ief) = obj.Fe(ief) - D(ief, jef) * eDOF_d(jef)*obj.system.sym_eps;          % eps?
                end
            end

            obj.F = Body.removeEps(obj.F, obj.system.sym_eps, obj.system.sym_eps_rot, true);
            obj.M = Body.removeEps(obj.M, obj.system.sym_eps, obj.system.sym_eps_rot, true);
            obj.Fe = Body.removeEps(obj.Fe, obj.system.sym_eps, obj.system.sym_eps_rot, true);
        end

        function calcGenForce(obj)
            % Constraint Forces can only be calculated once the whole
            % system loads were calculated. Therefore this is not possible
            % in prepareForces
            obj.applyContrLoads();

            obj.Fgen = - obj.v_p.' * (obj.F + obj.F_ext);
            obj.Fgen = obj.Fgen - obj.omega_p.' * (obj.M + obj.M_ext);            
            obj.Fgen = obj.Fgen - obj.e_p.' * (obj.Fe + obj.Fe_ext);
            obj.Fgen = simplify(obj.Fgen);
        end

        function applyContrLoads(obj)
            for i = 1:length(obj.children)
                if ~isempty(obj.ConstrLoads)
                    constr_forces = obj.system.getConstraintForce(obj.ConstrLoads{i}, false);
                    for j = 1:length(obj.ElasticDOF)
                        % M0 is already considered by the coupling of bodies via T_elast
                        e_force = obj.sid.frame(obj.ChildFrame(i)).phi.M1{j}' *constr_forces;
                        obj.Fe = obj.Fe + e_force*obj.ElasticDOF(j)*obj.system.sym_eps;
                    end
                end
            end
            obj.Fe = Body.removeEps(obj.Fe, obj.system.sym_eps, obj.system.sym_eps_rot, true);
        end
    end
end
