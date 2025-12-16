classdef ElasticFrame  < handle
    properties
        name (1,:) char = ''
        rframe (1,:) char = ''
        origin ElasticTaylor
        ap ElasticTaylor
        phi ElasticTaylor
        psi ElasticTaylor
        sigma ElasticTaylor
    end

    methods
        % Constructor
        function obj = ElasticFrame(nedof, name, rel_tol, system)
            arguments
                nedof
                name (1,:) char = ''
                rel_tol (1,1) double = 0
                system = []
            end
            if isstruct(nedof)
                frame_struct = nedof;
                
                obj.name = frame_struct.node;
                obj.rframe = frame_struct.rframe;
                obj.origin = ElasticTaylor(frame_struct.origin, 0, rel_tol, [name '_origin'], system);
                obj.ap = ElasticTaylor(frame_struct.AP, 0, rel_tol, [name '_ap'], system);
                obj.phi = ElasticTaylor(frame_struct.Phi, 0, rel_tol, [name '_phi'], system);
                obj.psi = ElasticTaylor(frame_struct.Psi, 0, rel_tol, [name '_psi'], system);
                obj.sigma = ElasticTaylor(frame_struct.sigma, 0, rel_tol, [name '_sigma'], system);                
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
    end
end
