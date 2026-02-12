classdef ElasticFrame  < handle
    properties
        name (1,:) char = ''
        rframe (1,:) char = ''
        origin ElasticTaylor    % M0: 3 x 1,  M1: 3 x 1 x nq;    M1=phi.M0 (table 6.9)
        ap ElasticTaylor        % M0: 3 x 3,  M1: 3 x 3 x nq;    M1: unit matrix, M1 = skew psi.M0 (table 6.9)
        phi ElasticTaylor       % M0: 3 x nq, M1: 3 x nq x nq
        psi ElasticTaylor       % M0: 3 x nq, M1: 3 x nq x nq
        sigma ElasticTaylor     % M0: 6 x 1,  M1: 6 x 1 x nq
        % according to standard, M1 dimensions 2 and 3 should be switched
    end

    methods
        % Constructor - supports two calling conventions
        %
        % Calling convention 1: Direct initialization
        %   obj = ElasticFrame(nedof, name)
        %
        %   nedof: number of elastic degrees of freedom
        %   name: name string for the frame (optional, default = '')
        %
        % Calling convention 2: Initialize from frame struct
        %   obj = ElasticFrame(frame_struct, name, tol, system)
        %
        %   frame_struct: struct with frame data (node, rframe, Phi, Psi, 
        %                 AP, sigma, origin)
        %   name: base name for generated parameters
        %   tol: relative tolerance for parameter generation (scalar or 
        %        [rel_tol, abs_tol])
        %   system: object with addParameter method callback
        %
        function obj = ElasticFrame(varargin)
            if nargin == 0
                error('ElasticFrame constructor requires at least one argument');
            end
            
            arg1 = varargin{1};
            
            if isstruct(arg1)
                % Calling convention 2: Initialize from struct
                frame_struct = arg1;
                name = '';
                tol = 0;
                system = [];
                
                if nargin >= 2
                    name = varargin{2};
                end
                if nargin >= 3
                    tol = varargin{3};
                end
                if nargin >= 4
                    system = varargin{4};
                end
                
                obj.name = frame_struct.node;
                obj.rframe = frame_struct.rframe;
                obj.phi = ElasticTaylor(frame_struct.Phi, 0, tol, [name '_phi'], system);

                % Compute origin.M1 from phi to reuse parameter symbols
                frame_struct.origin.M1 = arrayfun(@(i)obj.phi.M0(:, i), 1:size(obj.phi.M0, 2), UniformOutput=false);
                obj.origin = ElasticTaylor(frame_struct.origin, 0, tol, [name '_origin'], system);
                obj.psi = ElasticTaylor(frame_struct.Psi, 0, tol, [name '_psi'], system);

                % Compute AP.M1 from psi using cross-product matrix to reuse parameter symbols
                frame_struct.AP.M1 = arrayfun(@(i)crossmat(obj.psi.M0(:, i)), 1:size(obj.psi.M0, 2), UniformOutput=false);
                % Safeguard: frame node orientation must be identity matrix
                if ~all(frame_struct.AP.M0==eye(3))
                    error('Frame node orientation must always be a unit matrix.')
                end
                obj.ap = ElasticTaylor(frame_struct.AP, 0, tol, [name '_ap'], system);
                obj.sigma = ElasticTaylor(frame_struct.sigma, 0, tol, [name '_sigma'], system);
            else
                % Calling convention 1: Direct initialization
                nelastq = arg1;
                name = '';
                
                if nargin >= 2
                    name = varargin{2};
                end
                
                obj.name = name;
                obj.origin = ElasticTaylor(1, 3, 1, 0, nelastq, 0, 3);      % 3 x 1, M1: 3 x 1 x nelastq
                obj.ap = ElasticTaylor(1, 3, 3, 0, nelastq, 0, 3);          % 3 x 3, M1: 3 x 3 x nelastq
                obj.phi = ElasticTaylor(1, 3, nelastq, 0, nelastq, 0, 3);   % 3 x nelastq, M1: 3 x nelastq x nelastq
                obj.psi = ElasticTaylor(0, 3, nelastq, 0, nelastq, 0, 3);   % 3 x nelastq
                obj.sigma = ElasticTaylor(1, 6, 1, 0, nelastq, 0, 3);       % 6 x 1, M1: 6 x 1 x nelastq
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
