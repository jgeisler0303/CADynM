function matlabTemplateEngine(out_file_name, template_name, temp_data)
% matlabTemplateEngine(out_file_name, template_name, temp_data)
% Reads a template file and processes embedded MATLAB code blocks.
%
% Code blocks use the syntax:
%       {{ MATLAB_CODE }}
%
% The code is executed in a workspace where fields of 'temp_data'
% are available as variables. Output printed to stdout (fprintf,
% disp, etc.) is captured and substituted into the output stream.
%
% The processed output is written to 'out_file_name'.

    %----------------------------------------------------------------------
    % Read template
    %----------------------------------------------------------------------
    if ~exist(template_name, 'file')
        mydir= fileparts(mfilename("fullpath"));
        template_name = fullfile(mydir, 'templates', template_name);
    end
    template = fileread(template_name);

    %----------------------------------------------------------------------
    % Prepare output file
    %----------------------------------------------------------------------
    fid_out = fopen(out_file_name, 'w');
    if fid_out < 0
        error('Cannot open output file: %s', out_file_name);
    end
    cleanupObj = onCleanup(@() fclose(fid_out));

    % Write creation comment
    [~, ~, ext] = fileparts(out_file_name);
    [~, template_file_name, t_ext] = fileparts(template_name);
    switch ext
        case '.m'
            comment_char = '%%';
        case {'.c' '.cpp' '.h' '.hpp'}
            comment_char = '//';
    end

    d = datetime("today");
    fprintf(fid_out, '%s File created by matlabTemplateEngine on %s from %s\n\n', comment_char, string(d), [template_file_name t_ext]);

    %----------------------------------------------------------------------
    % Find all {{ ... }} blocks
    %----------------------------------------------------------------------
    %
    % Using regexp with tokens and positions:
    %   - 'start' = positions of '{{'
    %   - 'end'   = positions of '}}'
    %   - 'tokens' contains code inside
    %
    pattern = '\{\{(.*?)\}\}';  % non-greedy match
    [tokens, tokenStarts, tokenEnds] = regexp(template, pattern, ...
                                              'tokens', 'start', 'end');

    % If no code blocks, output entire template
    if isempty(tokens)
        fwrite(fid_out, template);
        return;
    end

    %----------------------------------------------------------------------
    % Process template by iterating through literal segments + code segments
    %----------------------------------------------------------------------
    currentPos = 1;

    for k = 1:numel(tokens)
        % Write literal text before this code block
        literalSegment = template(currentPos : tokenStarts(k)-1);
        fwrite(fid_out, literalSegment);

        % Extract code inside {{ ... }}
        code = strtrim(tokens{k}{1});

        % Execute code in a function workspace with temp_data unpacked
        outStr = executeCodeBlock(code, temp_data);

        % Write the captured output
        fwrite(fid_out, outStr);

        % Move current position
        currentPos = tokenEnds(k) + 1;
    end

    % Write final literal segment after last code block
    if currentPos <= length(template)
        fwrite(fid_out, template(currentPos:end));
    end
end


% =========================================================================
% Helper: Execute MATLAB code block with temp_data unpacked
% =========================================================================
function outStr = executeCodeBlock(code, temp_data)
    % Import values from struct into workspace of this function
    if isstruct(temp_data)
        fn = fieldnames(temp_data);
        for i = 1:numel(fn)
            assignin('caller', fn{i}, temp_data.(fn{i}));
        end
    end

    % Execute code and capture stdout (evalc captures printed output)
    try
        outStr = evalc(code);
    catch ME
        error('Error executing code block: %s\n%s', code, ME.message);
    end
end
