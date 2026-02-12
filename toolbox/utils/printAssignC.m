function ij_ = printAssignC(indent, name, var, type, selector)
arguments
    indent
    name {mustBeTextScalar}
    var 
    type {mustBeTextScalar}
    selector function_handle = @(x) true
end

[m, n] = size(var);
if ~ischar(indent)
    ind= repmat(' ', 1, indent);
else
    ind= indent;
end

ij= 0;
for i = 1:m
    for j = 1:n
        if ~selector(var(i, j)), continue, end
        switch type
            case 'numbered'
                index = sprintf('%d', i+(j-1)*m);
            case 'scalar'
                index = '';
            case 'vector'
                index = sprintf('[%d]', i-1+(j-1)*m);
            case 'vector_par'
                index = sprintf('(%d)', i-1+(j-1)*m);
            case 'matrix'
                index = sprintf('(%d, %d)', i-1, j-1);
            otherwise
                error('Unknown indexing "%s"', type)
        end

        code_str = ccode(var(i, j));
        pos = strfind(code_str, '=');
        if ~isempty(pos)
            code_str = code_str(pos+1:end); % remove "t0 = "
        end
        code_str = strrep(code_str, '_IDX', '(');
        code_str = strrep(code_str, 'XDI_', ')');
        code_str = strrep(code_str, 'PSTRUCT_', 'param.');
        if ~endsWith(code_str, ';')
            fprintf('%s%s%s = %s;\n', ind, name, index, code_str)
        else
            fprintf('%s%s%s = %s\n', ind, name, index, code_str)
        end
        ij = ij+1;
    end
end

if nargout>0
    ij_ = ij;
end