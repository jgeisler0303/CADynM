classdef SID  < handle
    properties
        name (1,:) char = ''
        comment (1,:) char = ''
        mass (1,1) = 0
        nelastq (1,1) double = 0 % number of elastic DOF (nq)
        ielastq (1,:) cell = {}
        frame (1,:) ElasticFrame = ElasticFrame.empty
        I ElasticTaylor         % M0: 3 x 3,  M1: 3 x 3 x nq
        md ElasticTaylor        % M0: 3 x 1,  M1: 3 x 1 x nq;     M1=reshape(Ct.M0', 3, 1, [])
        Ct ElasticTaylor        % M0: nq x 3, M1: nq x 3 x nq
        Cr ElasticTaylor        % M0: nq x 3, M1: nq x 3 x nq
        Me ElasticTaylor        % M0: nq x nq
        Gr ElasticTaylor        % M0: 3 x nq*3 -> 3 x 3 x nq (parameter numbering: nq x 3 x 3)
        Oe ElasticTaylor        % M0: nq x 6, M1: nq x 6 x nq 
        Ge ElasticTaylor        % M0: nq x nq*3 -> nq x 3 x nq (parameter numbering: nq x nq x 3)
        Ke ElasticTaylor        % M0: nq x nq
        De ElasticTaylor        % M0: nq x nq
        % according to standard, M1 dimensions 2 and 3 should be switched
    end

    methods
        % Constructor - supports two calling conventions
        %
        % Calling convention 1: Direct initialization
        %   obj = SID(nedof, n_frames, name)
        %
        %   nedof: number of elastic degrees of freedom
        %   n_frames: number of frames (optional, default = 0)
        %   name: name string for the SID object (optional, default = '')
        %
        % Calling convention 2: Initialize from SID struct
        %   obj = SID(sid_struct, tolerance, name, system)
        %
        %   sid_struct: struct with SID data (refmod, frame, Ct, Cr, md, I, 
        %               Me, Gr, Ge, Oe, Ke, De)
        %   tolerance: relative tolerance for parameter generation (scalar or 
        %              [rel_tol, abs_tol]). If non-zero, values are converted 
        %              to parameters
        %   name: base name for generated parameters
        %   system: object with addParameter method callback
        %
        function obj = SID(varargin)
            if nargin == 0
                error('SID constructor requires at least one argument');
            end
            
            arg1 = varargin{1};
            
            if isstruct(arg1)
                % Calling convention 2: Initialize from struct
                if nargin < 4
                    error('When using struct initialization, provide all required arguments: SID(sid_struct, tolerance, name, system)');
                end
                
                sid_struct = arg1;
                tol = varargin{2};
                name = varargin{3};
                system = varargin{4};
                
                obj.name = name;
                obj.nelastq = sid_struct.refmod.nelastq; % nq
                obj.comment = sid_struct.comment;
                obj.ielastq = sid_struct.refmod.ielastq;
                
                if ~all(tol==0)
                    obj.mass = system.addParameter([name '_mass'], [], sid_struct.refmod.mass);
                else
                    obj.mass = sid_struct.refmod.mass;
                end
    
                for i = 1:length(sid_struct.frame)
                    obj.frame(i) = ElasticFrame(sid_struct.frame(i), sprintf('%s_frame_%d', name, i), tol, system);
                end
                
                obj.Ct = ElasticTaylor(sid_struct.Ct, 0, tol, [name '_Ct'], system);
                obj.Cr = ElasticTaylor(sid_struct.Cr, 0, tol, [name '_Cr'], system);
                
                % Compute md.M1 from Ct to reuse parameter symbols 
                M = obj.Ct.M0.';
                sid_struct.md.M1 = arrayfun(@(i)M(:, i), 1:size(M, 2), UniformOutput=false);
                obj.md = ElasticTaylor(sid_struct.md, 0, tol, [name '_md'], system);
                
                obj.I = ElasticTaylor(sid_struct.I, 0, tol, [name '_I'], system);
                obj.Me = ElasticTaylor(sid_struct.Me, 0, tol, [name '_Me'], system);
                obj.Gr = ElasticTaylor(sid_struct.Gr, obj.nelastq, tol, [name '_Gr'], system);
                obj.Ge = ElasticTaylor(sid_struct.Ge, obj.nelastq, tol, [name '_Ge'], system);
                obj.Oe = ElasticTaylor(sid_struct.Oe, 0, tol, [name '_Oe'], system);
                obj.Ke = ElasticTaylor(sid_struct.Ke, 0, tol, [name '_K'], system);
                obj.De = ElasticTaylor(sid_struct.De, 0, tol, [name '_D'], system);
            else
                % Calling convention 1: Direct initialization
                nelastq = arg1; % number of elastic DOF (nq)
                n_frames = 0;
                name = '';
                
                if nargin >= 2
                    n_frames = varargin{2};
                end
                if nargin >= 3
                    name = varargin{3};
                end
                
                obj.name = name;
                obj.nelastq = nelastq;
                
                for i = 1:nelastq
                    obj.ielastq{i} = sprintf('Noname DOF%02d', i);
                end
    
                if n_frames > 0
                    for i = 1:n_frames
                        obj.frame(i) = ElasticFrame(nelastq, sprintf('frame%02d', i));
                    end
                end
                
                obj.md = ElasticTaylor(1, 3, 1, 0, nelastq, 0, 3);              % M0: 3 x 1,  M1: 3 x 1 x nq
                obj.I = ElasticTaylor(1, 3, 3, 0, nelastq, 0, 2);               % M0: 3 x 3,  M1: 3 x 3 x nq
                obj.Ct = ElasticTaylor(1, nelastq, 3, 0, nelastq, 0, 3);        % M0: nq x 3, M1: nq x 3 x nq
                obj.Cr = ElasticTaylor(1, nelastq, 3, 0, nelastq, 0, 3);        % M0: nq x 3, M1: nq x 3 x nq
                obj.Me = ElasticTaylor(0, nelastq, nelastq, 0, nelastq, 0, 2);  % M0: nq x nq
                obj.Gr = ElasticTaylor(0, 3, 3, nelastq, nelastq, 0, 3);        % M0: 3 x 3 x nq
                obj.Ge = ElasticTaylor(0, nelastq, 3, nelastq, nelastq, 0, 3);  % M0: nq x 3 x nq
                obj.Oe = ElasticTaylor(1, nelastq, 6, 0, nelastq, 0, 3);        % M0: nq x 6, M1: nq x 6 x nq
                obj.Ke = ElasticTaylor(0, nelastq, nelastq, 0, nelastq, 0, 2);  % M0: nq x nq
                obj.De = ElasticTaylor(0, nelastq, nelastq, 0, nelastq, 0, 2);  % M0: nq x nq
            end
        end

        % this function is only to compare the behavior of this new class
        % to the maxima implementation
        function write_maxima(obj, fn)
            fid = fopen(fn, 'w');
            fprintf(fid, '//* %d nodes, %d modes, generated by FEMBeam2SID MATLAB script on 16-Dec-2025 19:20:56 by jgeisler *//\n', length(obj.frame), obj.nelastq);
            fprintf(fid, '%s: emptyElasticMode(%d);\n\n', obj.name, obj.nelastq);

            fprintf(fid, '%s@refmod@mass: %s;\n', obj.name, char(obj.mass));
            fprintf(fid, 'blade@refmod@ielastq[1]: "Eigen Mode    1 :      0.693356 Hz";\n\n');

            for i = 1:length(obj.frame)
                obj.frame(i).write_maxima(fid, sprintf('frame[%d]', i), i)
            end

            obj.md.write_maxima(fid, [obj.name '@md']);
            obj.I.write_maxima(fid, [obj.name '@I']);
            obj.Ct.write_maxima(fid, [obj.name '@Ct']);
            obj.Cr.write_maxima(fid, [obj.name '@Cr']);
            obj.Me.write_maxima(fid, [obj.name '@Me']);
            obj.Gr.write_maxima(fid, [obj.name '@Gr']);
            obj.Ge.write_maxima(fid, [obj.name '@Ge']);
            obj.Oe.write_maxima(fid, [obj.name '@Oe']);
            obj.Ke.write_maxima(fid, [obj.name '@K']);
            obj.De.write_maxima(fid, [obj.name '@D']);

            fclose(fid);
        end
    end
end
