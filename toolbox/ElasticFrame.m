classdef ElasticFrame  < handle
    properties
        name (1,:) char = ''
        rframe (1,:) char = ''
        origin ElasticTaylor    % M0: 3x1,  M1: 3x1xnq;    M1=phi.M0 (table 6.9)
        ap ElasticTaylor        % M0: 3x3,  M1: 3x3xnq;    M1: unit matrix, M1 = skew psi.M0 (table 6.9)
        phi ElasticTaylor       % M0: 3xnq, M1: 3xnqxnq
        psi ElasticTaylor       % M0: 3xnq, M1: 3xnqxnq
        sigma ElasticTaylor     % M0: 6x1,  M1: 6x1xnq
        % according to standard, M1 dimensions 2 and 3 should be switched
    end

    methods
        % Constructor
        function obj = ElasticFrame(nedof, name, tol, system)
            arguments
                nedof
                name (1,:) char = ''
                tol double = 0
                system = []
            end
            if isstruct(nedof)
                frame_struct = nedof;
                
                obj.name = frame_struct.node;
                obj.rframe = frame_struct.rframe;
                obj.phi = ElasticTaylor(frame_struct.Phi, 0, tol, [name '_phi'], system);

                frame_struct.origin.M1 = reshape(obj.phi.M0, 3, 1, []);
                obj.origin = ElasticTaylor(frame_struct.origin, 0, tol, [name '_origin'], system);

                obj.psi = ElasticTaylor(frame_struct.Psi, 0, tol, [name '_psi'], system);

                frame_struct.AP.M1= sym.empty;
                for i= 1:frame_struct.AP.nq
                    frame_struct.AP.M1(:, :, i) = crossmat(obj.psi.M0(:, i));
                end
                % save-guard against formerly wrong sid computation in FEMBeam2SID
                if ~all(frame_struct.AP.M0==eye(3))
                    error('Frame node orientation must always be a unit matix.')
                end
                obj.ap = ElasticTaylor(frame_struct.AP, 0, tol, [name '_ap'], system);
                obj.sigma = ElasticTaylor(frame_struct.sigma, 0, tol, [name '_sigma'], system);                
            else
                if exist('name', 'var')
                    obj.name = name;
                end
    
                obj.origin = ElasticTaylor(1, 3, 1, 0, nedof, 0, 3);
                obj.ap = ElasticTaylor(1, 3, 3, 0, nedof, 0, 3);
                obj.phi = ElasticTaylor(1, 3, nedof, 0, nedof, 0, 3);
                obj.psi = ElasticTaylor(0, 3, nedof, 0, nedof, 0, 3);
                obj.sigma = ElasticTaylor(1, 6, 1, 0, nedof, 0, 3);
            end
        end

        % this function is only to compare the behavior of this new class
        % to the maxima implementation
        function write_maxima(obj, fid, name, i)
            fprintf(fid, '%s: emptyElasticFrame(%d, "node%d");\n', name, obj.origin.nq, i);

            obj.origin.write_maxima(fid, [name '@origin'])
            obj.ap.write_maxima(fid, [name '@ap'])
            obj.phi.write_maxima(fid, [name '@phi'])
            obj.psi.write_maxima(fid, [name '@psi'])
            obj.sigma.write_maxima(fid, [name '@sigma'])
        end
    end
end
