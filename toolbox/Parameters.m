classdef Parameters < handle

    properties (Access = private)
        param_syms struct = struct()
        param_dims struct = struct()
        param_struct struct = struct()
        param_values struct = struct()
    end

    methods
        function obj = Parameters(param_struct_)
            arguments
                param_struct_ struct = []
            end
            if ~isempty(param_struct_)
                obj.setParamRefStruct(param_struct_);
            end
        end

        function setParamRefStruct(obj, param_struct_)
            obj.param_struct = param_struct_;
        end

        function p = addParameter(obj, name, dims, value)
            arguments
                obj
                name
                dims = []
                value = []
            end
            if ~isempty(dims) && ~isempty(value) && ~all(dims==size(value))
                error('The dimensions parameter %s and the dimensions of the value %s don''t match.', mat2str(dims), mat2str(size(value)));
            end
            if isempty(dims) && ~isempty(value) && ~isscalar(value)
                dims = size(value);
            end
            if isfield(obj.param_struct, name)
                if ~isempty(value)
                    if ~all(size(obj.param_struct.(name))==size(value))
                        error('Dimension %s of parameter %s conflict with dimension of value in reference struct %s.', mat2str(size(value)), name, mat2str(size(obj.param_struct.(name))))
                    end
                else
                    value = obj.param_struct.(name);
                end
            end
            if isfield(obj.param_syms, name)
                % For now param symbols are always scalars. Vector or
                % matrix params are meant only as external parameters
                if ~all(obj.param_dims.(name)==size(value))
                    error('Trying to set value of dimension %s to parameter of dimension %s.', mat2str(size(value)), mat2str(obj.param_dims.(name)))
                end
            else
                if ~isempty(dims)
                    obj.param_syms.(name) = sym(name, 'real');
                    obj.param_dims.(name) = dims;
                else
                    obj.param_syms.(name) = sym(name, 'real');
                    obj.param_dims.(name) = [1 1];
                end
            end
            if ~isempty(value)
                obj.param_values.(name) = value;
            end
            if nargout>0
                p = obj.param_syms.(name);
            end
        end

        % Intercept property access
        function val = subsref(obj, S)
            if isscalar(S) && strcmp(S(1).type, '.')
                name = S(1).subs;

                % If field exists, return it
                if isfield(obj.param_syms, name)
                    val = obj.param_syms.(name);
                else
                    val = obj.addParameter(name, [], []);
                end
            else
                % Delegate all other indexing to builtin behavior
                if nargout>0
                    val = builtin('subsref', obj, S);
                else
                    builtin('subsref', obj, S);
                end
            end
        end

        % writing to a parameter sets the numeric value
        function obj = subsasgn(obj, S, value)
            if isscalar(S) && strcmp(S(1).type, '.')
                name = S(1).subs;

                obj.addParameter(obj, name, [], value)
            else
                obj = builtin('subsasgn', obj, S, value);
            end
        end

        function removeUnused(obj, vars)
            fn = fieldnames(obj.param_syms);
            for i = 1:length(fn)
                if ~ismember(fn{i}, vars)
                    obj.param_syms = rmfield(obj.param_syms, fn{i});
                    if isfield(obj.param_values, fn{i})
                        obj.param_values = rmfield(obj.param_values, fn{i});
                    end
                end
            end
        end

        function n = getNumParams(obj)
            n = length(fieldnames(obj.param_syms));
        end
        function n = getParamNames(obj)
            n = fieldnames(obj.param_syms);
        end
        function tf = paramInUse(obj, name)
            tf = isfield(obj.param_syms, name);
        end
        function tf = hasValue(obj, name)
            tf = isfield(obj.param_values, name) && ~isempty(obj.param_values.(name));
        end
        function v = getValue(obj, name)
            v = obj.param_values.(name);
        end
        function d = getDims(obj, name)
            d = obj.param_dims.(name);
        end
        function s = getParamSyms(obj)
            s = obj.param_syms;
        end
        function s = getParamValues(obj)
            s = obj.param_values;
        end
    end
end
