function model = makeCADynM(model_path, target_path, param, skip_gen, files_to_generate)

if ~iscell(files_to_generate)
    files_to_generate= {files_to_generate};
end

[~, model_name]= fileparts(model_path);
if startsWith(model_name, 'model')
    model_name = model_name(6:end);
end
if ~exist('target_path', 'var')
    target_path= fullfile(pwd, model_name);
end

if ~exist(target_path, 'dir')
    mkdir(target_path);
end

% Check which file to generate
skip_gen_file= false(length(files_to_generate), 1);
if ~exist('skip_gen', 'var') || isempty(skip_gen)
    dd_src= dir(model_path);
    for i= 1:length(files_to_generate)
        code_file= get_target_name(model_name, files_to_generate{i});
        code_path= fullfile(target_path, code_file);
        if exist(code_path, 'file')
            dd_dst= dir(code_path);
            if dd_src.datenum<dd_dst.datenum
                fprintf('Skipping generation of "%s"\n', code_file)
                skip_gen_file(i)= true;
            end
        end
    end
    skip_gen= all(skip_gen_file);
end

if ~skip_gen
    [model_dir, model_file] = fileparts(model_path);
    old_dir = pwd;
    cleanupObj = onCleanup(@()cd(old_dir));
    cd(model_dir);
    modelFunc = str2func(model_file);
    model = modelFunc(param);
    [~] = model.getEOM;
    model.removeUnusedParameters();

    for i= 1:length(files_to_generate)
        [code_file, template_name]= get_target_name(model_name, files_to_generate{i});
        
        target_file = fullfile(target_path, code_file);
        try
            matlabTemplateEngine(target_file, template_name, model)
        catch ME
            if exist(target_file, 'file')
                [p, n, e] = fileparts(target_file);
                movefile(target_file, fullfile(p, [n '_error' e]))
            end
            rethrow(ME)
        end
    end
else
    fprintf('Skipping CADynM generation because model is already up to date or you wanted it so\n');
    model = [];
end


function [final_name, template_name] = get_target_name(model_name, generator)
if generator(1)=='_'
    final_name= [model_name generator];
    template_name = generator(2:end);
else
    final_name= generator;
    template_name = generator;
end
template_name = [template_name '.mte'];