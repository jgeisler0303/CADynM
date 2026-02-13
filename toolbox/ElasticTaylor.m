classdef ElasticTaylor  < handle
    properties
        order (1,1) double      % order: here only 0 or 1
        nrow (1,1) double       % number of rows (origin, phi, psi, ap, md, I, Gr: 3; sigma: 6; Ct, Cr, Me, Ge, Oe, K, D: nq)
        ncol (1,1) double       % number of columns (origin, sigma, md: 1; ap, I, Ct, Cr: 3; Oe: 6; phi, psi, Me, K, D: nq; Gr, Ge: nq*3)
        nelem (1,1) double      % number of elements 0 or nq for Gr and Ge
        nq (1,1) double         % number of elastic DOF
        nqn (1,1) double        % dimension of order >1, unused here
        structure (1,1) double  % 0: all zeros; 1: M is diagonal; 2: M is symmetric; 3: M0 is nrow x ncol x nelem; 4: not implemented
        M0  % nrow x ncol or nelem cell array of nrow x ncol matrices
        M1  % nq cell array of nrow x ncol matrices (according to specification table 6.2: nrow x nq x ncol)
    end

    methods
        % Constructor - supports two calling conventions
        %
        % Calling convention 1: Direct initialization
        %   obj = ElasticTaylor(order, nrow, ncol, nelem, nq, nqn, structure)
        %
        %   order: Taylor expansion order (0 or 1)
        %   nrow: number of rows
        %   ncol: number of columns
        %   nelem: number of elements
        %   nq: number of elastic DOF
        %   nqn: dimension of order > 1
        %   structure: structure type
        %             0: all zeros
        %             1: diagonal matrix
        %             2: symmetric matrix
        %             3: general matrix
        %             4: unit matrix
        % M0: nrow x ncol matrix or nrow x ncol x nelem if nelem>0
        % M1: cell array of nrow x ncol matrices, length nq
        %
        % Calling convention 2: Initialize from struct
        %   obj = ElasticTaylor(taylor_struct, nelem, tol, name, system)
        %
        %   taylor_struct: struct with fields (order, nrow, ncol, nq, nqn, 
        %                  structure, M0, M1)
        %   nelem: number of elements
        %   tol: tolerance (scalar relative tolerance or [rel_tol, abs_tol])
        %   name: base name for generated parameters
        %   system: system object with addParameter method callback
        %
        function obj = ElasticTaylor(varargin)
            if nargin < 3
                error('At least three arguements required for construction.');
            elseif isstruct(varargin{1})
                % Calling convention 2: Initialize from struct
                if nargin~=5
                    error('When initializing from struct, five arguments are needed: taylor_struct, nelem, tol, name, system')
                end
                taylor_struct = varargin{1};
                nelem = varargin{2};
                tol = varargin{3};
                name = varargin{4};
                system = varargin{5};
                
                obj.order = taylor_struct.order;
                obj.nrow = taylor_struct.nrow;
                obj.ncol = taylor_struct.ncol;
                obj.nq = taylor_struct.nq;
                obj.nqn = taylor_struct.nqn;
                obj.structure = taylor_struct.structure;
                obj.nelem = nelem;

                if isempty(tol)
                    rel_tol = 0;
                    abs_tol = 0;
                elseif isscalar(tol)
                    rel_tol = tol;
                    abs_tol = inf;
                else
                    rel_tol = tol(1);
                    abs_tol = tol(2);
                end
                tol_ = min(abs_tol, norm(taylor_struct.M0(:))*abs(rel_tol));

                switch obj.structure
                    case 0
                        if ~all(abs(taylor_struct.M0)>abs_tol)
                            error('Inconsistent structure 0 with non-zero M0 in ElasticTaylor initialization.')
                        else
                            taylor_struct.M0 = taylor_struct.M0*0;
                        end
                    case 1
                        M0_ref = diag(diag(taylor_struct.M0));
                        if any(abs(taylor_struct.M0-M0_ref)>tol_)
                            error('Inconsistent structure 1 with non-diagonal M0 in ElasticTaylor initialization.')
                        else
                            taylor_struct.M0 = M0_ref;
                        end
                    case 2
                        if any(abs(taylor_struct.M0.' - taylor_struct.M0)>tol_)
                            error('Inconsistent structure 2 with non-symmetric M0 in ElasticTaylor initialization.')
                        else
                            taylor_struct.M0 = 0.5*(taylor_struct.M0 + taylor_struct.M0.');
                        end
                    case 4
                        M0_ref = eye(size(taylor_struct.M0));
                        if any(abs(taylor_struct.M0.' - taylor_struct.M0)>tol_)
                            error('Inconsistent structure 4 with non-unit M0 in ElasticTaylor initialization.')
                        else
                            taylor_struct.M0 = M0_ref;
                        end
                end
                if nelem>0 && ~iscell(taylor_struct.M0)
                    % if M0 is given as a 3D array, reshape it to a cell array of 2D matrices
                    % if it already is a cell array, we assume it is correctly formatted as nrow x ncol x nelem
                    M0_ = reshape(taylor_struct.M0, size(taylor_struct.M0, 1), size(taylor_struct.M0, 2)/nelem, nelem);
                else
                    M0_ = taylor_struct.M0;
                end
                obj.M0 = obj.setElements(M0_, 0, tol, [name '0'], system);

                % TODO: check structure for M1
                if obj.order>0
                    obj.M1 = obj.setElements(taylor_struct.M1, 1, tol, [name '1'], system);
                end
            else
                % Calling convention 1: Direct initialization with individual parameters
                if nargin~=5
                    error('When initializing empty ElasticTaylor, seven arguments are needed: order, nrow, ncol, nelem, nq, nqn, structure.')
                end
                obj.order = varargin{1};
                obj.nrow = varargin{2};
                obj.ncol = varargin{3};
                obj.nelem = varargin{4};
                obj.nq = varargin{5};
                obj.nqn = varargin{6};
                obj.structure = varargin{7};
                
                if obj.nelem>1
                    obj.M0 = repmat({zeros(obj.nrow, obj.ncol)}, obj.nelem);
                else
                    obj.M0 = zeros(obj.nrow, obj.ncol);
                end
                
                if obj.order>0
                    if obj.nelem>0
                        error('Order 1 with nelem>0 not implemented yet.')
                    else
                        obj.M1 = repmat({zeros(obj.nrow, obj.ncol)}, obj.nq);
                    end                        
                end
                        
                if obj.order>1
                    error("Elements of class taylor with order > 1 currently not supported")
                end
            end
        end

        % evaluate the taylor data at a given point
        function r = evalTaylor(obj, dof, eps)
            r = obj.M0;
            % TODO: correct dimenations?
            if obj.order>0
                for i = 1:obj.nq
                    if obj.nelem>0
                        % r = r + obj.M1(:, :, :, i)*dof(i) * eps; % eps correct to use here?
                        error('not finished ElasticTaylor access to M1 with 4 dimensions.')
                    else
                        r = r + obj.M1{i}*dof(i) * eps;
                    end
                end
            end
        end

        % this function is only to compare the behavior of this new class
        % to the maxima implementation
        function write_maxima(obj, fid, name)
            if obj.nelem==0
                for i = 1:size(obj.M0, 1)
                    for j = 1:size(obj.M0, 2)
                        if isempty(symvar(obj.M0(i, j))), continue, end
                        fprintf(fid, '%s@M0[%d, %d]: %s;\n', name, i, j, char(obj.M0(i, j)));
                    end
                end
            else
                for i = 1:size(obj.M0, 1)
                    for j = 1:size(obj.M0, 2)
                        for k = 1:size(obj.M0, 3)
                            if isempty(symvar(obj.M0(i, j, k))), continue, end
                            fprintf(fid, '%s@M0[%d][%d, %d]: %s;\n', name, k, i, j, char(obj.M0(i, j, k)));
                        end
                    end
                end
            end
            fprintf(fid, '\n');
            if obj.order>0
                for k = 1:size(obj.M1, 3)
                    for i = 1:size(obj.M1, 1)
                        for j = 1:size(obj.M1, 2)
                            if isempty(symvar(obj.M1(i, j, k))), continue, end
                            fprintf(fid, '%s@M1[%d][%d, %d]: %s;\n', name, k, i, j, char(obj.M1(i, j, k)));
                        end
                    end
                end            
            end
            fprintf(fid, '\n');
        end
    end

    methods (Access = private)
        function M = setElements(obj, M_, order, tol, name, system)
            if iscell(M_) || isa(M_, 'sym') || isa(M_, 'msym')
                % if M_ is already a cell array or symbolic, we assume it was already processed and just return it as M
                M = M_;
                return
            end

            if isscalar(tol)
                rel_tol = tol;
                abs_tol = inf;
            else
                rel_tol = tol(1);
                abs_tol = tol(2);
            end

            % rel_tol==0 means we keep the numbers
            if rel_tol==0
                if order==1 && obj.nelem>0
                    error('4-dimensional M1 not supported in ElasticTaylor initialization.')
                    % [R,C] = meshgrid(1:size(M_, 3), 1:size(M_, 4));
                    % M = arrayfun(@(r,c)M_(:,:,r,c), R, C, UniformOutput=false);
                elseif order==1 || obj.nelem>0
                    M = arrayfun(@(r)M_(:,:,r), 1:size(M_, 3), UniformOutput=false);
                else
                    % M0 with nelem=0
                    M = M_;
                end
            else
            % all numbers > tol are turned into parameters, the rest is set to zero                
                tol = min(abs_tol, norm(M_(:))*abs(rel_tol));
                % identify the entries that should be turned into parameters (params_idx is a logical index)
                if rel_tol<0
                    % rel_tol<0 means numbers that are exactly 1 or -1 are not turned into parameters
                    params_idx= abs(M_)>tol & abs(M_)~=1;
                else
                    params_idx= abs(M_)>tol;
                end
                M_(abs(M_)<=tol) = 0;

                % prepare the symbolic M for parameter replacement
                if order==1 && obj.nelem>0
                    error('4-dimensional M1 not supported in ElasticTaylor initialization.')
                    % [R,C] = meshgrid(1:size(M_, 3), 1:size(M_, 4));
                    % M = arrayfun(@(r,c)system.sym(M_(:,:,r,c)), R, C, UniformOutput=false);
                elseif order==1 || obj.nelem>0
                    % M0 with nelem>0 or M1 with nelem=0
                    M = arrayfun(@(r)system.sym(M_(:,:,r)), 1:size(M_, 3), UniformOutput=false);
                else
                    % M0 with nelem=0
                    M = system.sym(M_);
                end

                params_idx_lin = find(params_idx);
                dims = obj.checkDimensions(M_, order);

                % convert linear indices to subscripts for parameter naming
                [subs{1:length(dims)}] = ind2sub(dims, params_idx_lin(:));
                subsMat = [subs{:}];

                % replace the identified entries with parameters
                for i = 1:length(params_idx_lin)
                    if obj.structure==2
                        % reuses parameters of symmetric entries
                        if subsMat(i, 1)<subsMat(i, 2)
                            swap = subsMat(i, 2);
                            subsMat(i, 2) = subsMat(i, 1);
                            subsMat(i, 1) = swap;
                        end
                    end
                    if order==0 && obj.nelem==0
                        pname = [name sprintf('_%d', subsMat(i, :))];
                    else
                        pname = [name sprintf('_%d', [subsMat(i, end) subsMat(i, 1:end-1)])];
                    end
                    p = system.addParameter(pname, [], M_(params_idx_lin(i)));
                    if iscell(M)
                        M{subsMat(i, 3)}(subsMat(i, 1), subsMat(i, 2)) = p;
                    else
                        M(subsMat(i, 1), subsMat(i, 2)) = p;
                    end
                end
            end
        end

        function dims = checkDimensions(obj, M, order)
            if order==0
                if obj.nelem~=0
                    dims = [obj.nrow, obj.ncol/obj.nelem, obj.nelem];
                else
                    dims = [obj.nrow, obj.ncol];
                end
            else
                if obj.nelem~=0
                    dims = [obj.nrow, obj.ncol/obj.nelem, obj.nelem, obj.nq];
                else
                    dims = [obj.nrow, obj.ncol, obj.nq];
                end
            end
            if ndims(M)>length(dims) || ~isequal(size(M, 1:length(dims)), dims)
                error('Size of M does not match expected dimensions based on nrow, ncol, nelem, and nq.')
            end
        end
    end
end
