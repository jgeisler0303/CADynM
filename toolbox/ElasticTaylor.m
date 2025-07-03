classdef ElasticTaylor  < handle
    properties
        order (1,1) double
        nrow (1,1) double
        ncol (1,1) double
        nelem (1,1) double
        nq (1,1) double
        nqn (1,1) double
        structure (1,1) double % 0: no M0; 1: M is diagonla; 2: M is symmetric; 3: M0 is nrow x ncol x nelem; 4: not implemented
        M0 double
        M1 double % row x ncol x nq
    end

    methods
        % Constructor
        function obj = ElasticTaylor(order, nrow, ncol, nelem, nq, nqn, structure)
            if isstruct(order)
                taylor_struct = order;
                if exist('nrow', 'var')
                    nelem = nrow;
                else
                    nelem= 0;
                end

                obj.order = taylor_struct.order;
                obj.nrow = taylor_struct.nrow;
                obj.ncol = taylor_struct.ncol;
                obj.nq = taylor_struct.nq;
                obj.nqn = taylor_struct.nqn;
                obj.structure = taylor_struct.structure;
                
                obj.nelem = nelem;
                if nelem>0
                    obj.M0 = reshape(taylor_struct.M0, size(taylor_struct.M0, 1), size(taylor_struct.M0, 2)/nelem, nelem);
                else
                    obj.M0 = taylor_struct.M0;
                end

                if isfield('M1', taylor_struct)
                    obj.M1 = taylor_struct.M1;
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
end
