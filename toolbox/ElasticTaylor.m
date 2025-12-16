classdef ElasticTaylor  < handle
    properties
        order (1,1) double
        nrow (1,1) double
        ncol (1,1) double
        nelem (1,1) double
        nq (1,1) double
        nqn (1,1) double
        structure (1,1) double % 0: no M0; 1: M is diagonla; 2: M is symmetric; 3: M0 is nrow x ncol x nelem; 4: not implemented
        M0 
        M1  % row x ncol x nq
    end

    methods
        % Constructor
        function obj = ElasticTaylor(order, nrow, ncol, nelem, nq, nqn, structure)
            arguments
                order
                nrow (1,1) double = 0
                ncol (1,1) double = 0 
                nelem = []
                nq = []
                nqn = []
                structure = []
            end
            if isstruct(order)
                taylor_struct = order;
                rel_tol = ncol;
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

                % TODO: fix dimensions
                if nelem>0
                    M0_ = reshape(taylor_struct.M0, size(taylor_struct.M0, 1), size(taylor_struct.M0, 2)/nelem, nelem);
                else
                    M0_ = taylor_struct.M0;
                end
                obj.M0 = obj.setElements(M0_, rel_tol, [name '0'], system);

                if isfield(taylor_struct, 'M1')
                    obj.M1 = obj.setElements(taylor_struct.M1, rel_tol, [name '1'], system);
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
    end
    methods (Static, Access = private)
        function M = setElements(M_, rel_tol, name, system)
            if rel_tol==0
                M = M_;
            else
                tol = norm(M_)*abs(rel_tol);
                if rel_tol<0
                    params_idx= abs(M_)>tol & abs(M_)~=1;
                else
                    params_idx= abs(M_)>tol;
                end
                M_(abs(M_)<=tol) = 0;
                M = sym(M_);

                params_idx_lin = find(params_idx);
                if isvector(M_)
                    subsMat = params_idx_lin;
                else
                    [subs{1:ndims(M_)}] = ind2sub(size(M_), params_idx_lin);
                    subsMat = [subs{:}];
                end
                for i = 1:size(subsMat, 1)
                    pname = [name sprintf('_%d', subsMat(i, :))];
                    p = system.addParameter(pname);
                    M(params_idx_lin(i)) = p;
                end
            end
        end
    end
end
