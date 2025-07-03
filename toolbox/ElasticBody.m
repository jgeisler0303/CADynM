classdef ElasticBody  < Body
    properties
        ElasticDOF (:,1) sym = []
        sid SID
        Fe (:, 1) sym = 0
        Fe_ext (:, 1) sym = 0
        e_p (:,:) sym = []
    end

    methods
        % Constructor
        function obj = ElasticBody(eDOF, sid, name, description)
            if ~isempty(eDOF) && length(eDOF)~=sid.refmod.nelastq
                error('Number of elastic DOF (%d) and DOF in SID (%d) must match.', length(eDOF), sid.nelastq)
            end
            for i = 1:length(eDOF)
                obj.ElasticDOF(i) = eDOF;
            end
            
            obj.sid = sid;
            if nargin > 2
                obj.Name = name;
            end
            if nargin > 3
                obj.Description = description;
            end
        end

        function prepareKinematics(obj)
            obj.e_p= jacobian(obj.ElasticDOF, obj.system.q);
            prepareKinematicsBase(obj)
        end

        function prepareForces(obj)
            sym_eps = sym('eps', 'real');
            eDOF_d= diff(obj.ElasticDOF, obj.system.time);
            eDOF_dd= diff(eDOF_d, obj.system.time);
            
            rotationToGlobal = obj.T0(1:3,1:3);
            rotationToLocal = rotationToGlobal.';
            
            obj.F = obj.sid.mass * (obj.system.gravity - obj.a0);
            
            md_global = rotationToGlobal * obj.sid.md.evalTaylor(obj.ElasticDOF, sym_eps);
            obj.F = obj.F - crossmat(obj.alpha0) * md_global;
            obj.F = obj.F - crossmat(obj.omega0) * (crossmat(obj.omega0) * md_global);

            Ct = evalTaylor(obj.sid.Ct, obj.ElasticDOF, sym_eps);
            Ct_d = zeros(3, 1);
            Ct_dd = zeros(3, 1);
            for ief = 1:length(obj.ElasticDOF)
                Ct_d = Ct_d + Ct(ief, :).' * eDOF_d(ief)*sym_eps;
                Ct_dd = Ct_dd + Ct(ief, :).' * eDOF_dd(ief)*sym_eps;
            end
            obj.F = obj.F - 2 * crossmat(obj.omega0) * rotationToGlobal * Ct_d;
            obj.F = obj.F - rotationToGlobal * Ct_dd;

            Phi = evalTaylor(obj.sid.I, obj.ElasticDOF, sym_eps);

            PhiG_global = rotationToGlobal * Phi * rotationToLocal;

            obj.M = -PhiG_global * obj.alpha0 - crossmat(obj.omega0) * (PhiG_global * obj.omega0);

            omega_local = rotationToLocal * obj.omega0;

            obj.M = obj.M + crossmat(md_global) * (obj.system.gravity - obj.a0);

            Cr = evalTaylor(obj.sid.Cr, obj.ElasticDOF, sym_eps);
            Cr_dd = zeros(3, 1);
            for ief = 1:length(obj.ElasticDOF)
                Cr_dd = Cr_dd + Cr(ief, :).' * eDOF_dd(ief)*sym_eps;
            end
            obj.M = obj.M - rotationToGlobal * Cr_dd;

            Gr = evalTaylor(obj.sid.Gr, obj.ElasticDOF, sym_eps);
            Gr_ = zeros(3, 3);
            for ief = 1:length(obj.ElasticDOF) 
                Gr_ = Gr_ + Gr(:, :, ief) * eDOF_d(ief)*sym_eps;
            end
            obj.M = obj.M - rotationToGlobal * Gr_ * omega_local;

            grav_local = rotationToLocal * (obj.a0 - obj.system.gravity);
            alpha_local = rotationToLocal * obj.alpha0;

            obj.Fe = zeros(length(obj.ElasticDOF), 1);

            for ief = 1:length(obj.ElasticDOF)
                obj.Fe(ief) = - Ct(ief, :) * grav_local;
                obj.Fe(ief) = obj.Fe(ief) - Cr(ief, :) * alpha_local;
                for jef = 1:length(obj.ElasticDOF)
                    obj.Fe(ief) = obj.Fe(ief) - obj.sid.Me.M0(ief, jef) * eDOF_dd(jef)*sym_eps;
                end
            end

            w = [omega_local(1)*omega_local(1), omega_local(2)*omega_local(2), omega_local(3)*omega_local(3), omega_local(1)*omega_local(2), omega_local(2)*omega_local(3), omega_local(1)*omega_local(3)];
            obj.Fe = obj.Fe - evalTaylor(obj.sid.Oe, obj.ElasticDOF, sym_eps) * w.';

            Ge = evalTaylor(obj.sid.Ge, obj.ElasticDOF, sym_eps);
            for ief = 1:length(obj.ElasticDOF)
                Ge_ = zeros(1, 3);
                for jef = 1:length(obj.ElasticDOF)
                    % TODO = check for correct order of jef and ief
                    Ge_ = Ge_ + Ge(ief, :, jef) * eDOF_d(jef)*sym_eps;
                end
                obj.Fe(ief) = obj.Fe(ief) - Ge_ * omega_local;
            end

            K = evalTaylor(obj.sid.Ke, obj.ElasticDOF, sym_eps);
            D = evalTaylor(obj.sid.De, obj.ElasticDOF, sym_eps);
            for ief = 1:length(obj.ElasticDOF)
                for jef = 1 : length(obj.ElasticDOF)
                    obj.Fe(ief) = obj.Fe(ief) - K(ief, jef) * obj.ElasticDOF(jef)*sym_eps;  % eps?
                    obj.Fe(ief) = obj.Fe(ief) - D(ief, jef) * eDOF_d(jef)*sym_eps;          % eps?
                end
            end
        end

        function calcGenForce(obj)
            obj.Fgen = obj.v_p.' * (obj.F + obj.F_ext);
            obj.Fgen = obj.Fgen + obj.omega_p.' * (obj.M + obj.M_ext);            
            obj.Fgen = obj.Fgen + obj.e_p.' * (obj.Fe + obj.Fe_ext);           
        end
    end
end
% Trot_elast(sys, nbody, nframe):= block(
%     ident(3) + apply("+", makelist(eps*eps_rot*sys@sid[nbody]@frame[nframe]@ap@M1[i]*sys@states[sys@elastic_dof_idx[nbody][i]], i, 1, length(sys@elastic_dof_idx[nbody])))
% );
% 
% Telast(nbody, nframe):= block([Toi, xyz, ori, z_elast],
%     z_elast: [concat(z_elast_, nbody, '_, nframe)[1], concat(z_elast_, nbody, '_, nframe)[2], concat(z_elast_, nbody, '_, nframe)[3]],
%     xyz: sid[nbody]@frame[nframe]@origin@M0 + apply("+", makelist(eps*sid[nbody]@frame[nframe]@origin@M1[i]*elastic_dof[nbody][i], i, 1, length(elastic_dof[nbody]))),
%     ori: ident(3) + apply("+", makelist(eps*eps_rot*sid[nbody]@frame[nframe]@ap@M1[i]*elastic_dof[nbody][i], i, 1, length(elastic_dof[nbody]))),
% 
%     z_list: append(z_list, z_elast),
%     node_forces: endcons([nbody, nframe, concat(z_elast_, nbody, '_, nframe)], node_forces), 
%     Toi: Tdisp(z_elast[1], z_elast[2], z_elast[3]) . Tdisp(xyz[1, 1], xyz[2, 1], xyz[3, 1]) . addcol(addrow(ori, [0, 0, 0]), [0, 0, 0, 1])
% );
