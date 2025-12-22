classdef ElasticTaylor  < handle
    properties
        order (1,1) double      % order: here only 0 or 1
        nrow (1,1) double       % number of rows (origin, phi, psi, ap, md, I, Gr: 3; sigma: 6; Ct, Cr, Me, Ge, Oe, K, D: nq)
        ncol (1,1) double       % number of columns (origin, sigma, md: 1; ap, I, Ct, Cr: 3; Oe: 6; phi, psi, Me, K, D: nq; Gr, Ge: nq*3)
        nelem (1,1) double      % number of elements 0 or nq for Gr and Ge
        nq (1,1) double         % number of elastic DOF
        nqn (1,1) double        % dimension of order >1
        structure (1,1) double  % 0: all zeros; 1: M is diagonal; 2: M is symmetric; 3: M0 is nrow x ncol x nelem; 4: not implemented
        M0  % nrow x ncol
        M1  % according to specification table 6.2: nrow x nq x ncol
    end

    methods
        % Constructor
        % TODO: call via varargin
        % alternative call:
        % ElasticTaylor(taylor_struct, nelem, tol, name, system)
        % tol can be relative tolerance or [rel_tol, abs_tol]
        % system allows call-back to create parameters
        function obj = ElasticTaylor(order, nrow, ncol, nelem, nq, nqn, structure)
            arguments
                order
                nrow (1,1) double = 0
                ncol double = 0 
                nelem = []
                nq = []
                nqn = []
                structure = []          % 0: zero matrix, 1: diagonal matrix, 2: symmetric matrix, 3: any matriy, 4: unit matrix
            end
            if isstruct(order)
                taylor_struct = order;
                tol = ncol;
                name = nelem;
                system = nq;
                nelem = nrow;

                obj.order = taylor_struct.order;
                obj.nrow = taylor_struct.nrow;
                obj.ncol = taylor_struct.ncol;
                obj.nq = taylor_struct.nq;
                obj.nqn = taylor_struct.nqn;
                obj.structure = taylor_struct.structure;
                
                obj.nelem = nelem;

                switch obj.structure
                    case 0
                        M0_ = taylor_struct.M0*0;
                    case 1
                        M0_ = diag(diag(taylor_struct.M0));
                    case 2
                        M0_ = 0.5*(taylor_struct.M0+taylor_struct.M0');
                    case 3
                        M0_ = taylor_struct.M0;
                    case 4
                        M0_ = eye(size(taylor_struct.M0));
                end
                if nelem>0
                    M0_ = reshape(M0_, size(taylor_struct.M0, 1), size(taylor_struct.M0, 2)/nelem, nelem);
                end
                obj.M0 = obj.setElements(M0_, 0, tol, [name '0'], system);

                % TODO: consider structure for M1
                if obj.order>0
                    if isa(taylor_struct.M1, 'sym')
                        obj.M1 = taylor_struct.M1;
                    else
                        obj.M1 = obj.setElements(taylor_struct.M1, 1, tol, [name '1'], system);
                    end
                end
            else
                obj.order = order;
                obj.nrow = nrow;
                obj.ncol = ncol;
                obj.nelem = nelem;
                obj.nq = nq;
                obj.nqn = nqn;
                obj.structure = structure;
                
                if nelem>1
                    obj.M0 = zeros(nrow, ncol, nelem);
                else
                    obj.M0 = zeros(nrow, ncol);
                end
                
                if order>0
                    if nelem>1
                        obj.M1 = repmat(obj.M0, 1, 1, 1, nq);
                    else
                        obj.M1 = repmat(obj.M0, 1, 1, nq);
                    end                        
                end
                        
                if order>1
                    error("Elements of class taylor with order > 1 currently not supported")
                end
            end
        end

        % evaluate the taylor data at a given point
        function r = evalTaylor(obj, dof, eps)
            r = obj.M0;
            % TODO: correct dimenations?
            if obj.order>0 && ~isempty(obj.M1)
                for i = 1:obj.nq
                    if obj.nelem>0
                        r = r + obj.M1(:, :, :, i)*dof(i) * eps; % eps correct to use here?
                    else
                        r = r + obj.M1(:, :, i)*dof(i) * eps; % eps correct to use here?
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
            if isscalar(tol)
                rel_tol = tol;
                abs_tol = inf;
            else
                rel_tol = tol(1);
                abs_tol = tol(2);
            end
            if rel_tol==0
                M = M_;
            else
                tol = min(abs_tol, norm(M_(:))*abs(rel_tol));
                if rel_tol<0
                    params_idx= abs(M_)>tol & abs(M_)~=1;
                else
                    params_idx= abs(M_)>tol;
                end
                M_(abs(M_)<=tol) = 0;
                M = sym(M_);

                params_idx_lin = find(params_idx);
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
                [subs{1:length(dims)}] = ind2sub(dims, params_idx_lin(:));
                subsMat = [subs{:}];

                for i = 1:length(params_idx_lin)
                    if obj.structure==2
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
                    M(params_idx_lin(i)) = p;
                end
            end
        end
    end
end
