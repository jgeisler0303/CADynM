classdef SID  < handle
    properties
        comment (1,:) char = ''
        mass (1,1) = 0
        nelastq (1,1) double = 0
        ielastq (1,:) cell = {}
        frame (1,:) ElasticFrame = ElasticFrame.empty
        I ElasticTaylor 
        md ElasticTaylor
        Ct ElasticTaylor
        Cr ElasticTaylor
        Me ElasticTaylor
        Gr ElasticTaylor
        Oe ElasticTaylor
        Ge ElasticTaylor
        Ke ElasticTaylor
        De ElasticTaylor
    end

    methods
        % Constructor
        function obj = SID(nedof, n_frames, name, system)
            arguments
                nedof
                n_frames (1,1) double = 0;
                name = ''
                system = []
            end
            if isstruct(nedof)
                sid_struct= nedof;
                % Specify relative tolerance for which parameters shall be
                % omitted.
                % If non-zero, values are changed to parameters.
                rel_tol = n_frames;

                obj.nelastq = sid_struct.refmod.nelastq;
                obj.comment = sid_struct.comment;
                obj.ielastq = sid_struct.refmod.ielastq;
                if rel_tol~=0
                    obj.mass= system.addParameter([name '_mass'], [], sid_struct.refmod.mass);
                else
                    obj.mass= sid_struct.refmod.mass;
                end
    
                for i = 1:length(sid_struct.frame)
                    obj.frame(i) = ElasticFrame(sid_struct.frame(i), sprintf('%s_frame_%d', name, i), rel_tol, system);
                end
                
                obj.md = ElasticTaylor(sid_struct.md, 0, rel_tol, [name '_md'], system);
                obj.I = ElasticTaylor(sid_struct.I, 0, rel_tol, [name '_I'], system); % 2nd order element currently not supported
                obj.Ct = ElasticTaylor(sid_struct.Ct, 0, rel_tol, [name '_Ct'], system);
                obj.Cr = ElasticTaylor(sid_struct.Cr, 0, rel_tol, [name '_Cr'], system);
                obj.Me = ElasticTaylor(sid_struct.Me, 0, rel_tol, [name '_Me'], system);
                obj.Gr = ElasticTaylor(sid_struct.Gr, obj.nelastq, rel_tol, [name '_Gr'], system);
                obj.Ge = ElasticTaylor(sid_struct.Ge, obj.nelastq, rel_tol, [name '_Ge'], system);
                obj.Oe = ElasticTaylor(sid_struct.Oe, 0, rel_tol, [name '_Oe'], system);
                obj.Ke = ElasticTaylor(sid_struct.Ke, 0, rel_tol, [name '_K'], system);
                obj.De = ElasticTaylor(sid_struct.De, 0, rel_tol, [name '_D'], system);
            else
                obj.nelastq = nedof;
                for i = 1:nedof
                    obj.ielastq{i} = sprintf('Noname DOF%02d', i);
                end
    
                
                if exist('n_frames', 'var')
                    for i = 1:n_frames
                        obj.frame(i) = ElasticFrame(nedof, sprintf('frame%02d', i));
                    end
                end
                
                obj.md = ElasticTaylor(1, 3, 1, 0, nedof, 0, 3);
                obj.I = ElasticTaylor(1, 3, 3, 0, nedof, 0, 2); % 2nd order element currently not supported
                obj.Ct = ElasticTaylor(1, nedof, 3, 0, nedof, 0, 3);
                obj.Cr = ElasticTaylor(1, nedof, 3, 0, nedof, 0, 3);
                obj.Me = ElasticTaylor(0, nedof, nedof, 0, nedof, 0, 2);
                obj.Gr = ElasticTaylor(0, 3, 3, nedof, nedof, 0, 3);
                obj.Ge = ElasticTaylor(0, nedof, 3, nedof, nedof, 0, 3);
                obj.Oe = ElasticTaylor(1, nedof, 6, 0, nedof, 0, 3);
                obj.Ke = ElasticTaylor(1, nedof, nedof, 0, nedof, 0, 2);
                obj.De = ElasticTaylor(0, nedof, nedof, 0, nedof, 0, 2);
            end
        end
    end
end
