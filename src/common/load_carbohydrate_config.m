function config = load_carbohydrate_config(repo_dir, molecule)
% Load a molecule-specific JSON workflow configuration.
% Paths are resolved relative to repo_dir; no user-machine paths are stored.
if nargin < 1 || isempty(repo_dir)
    this_file = mfilename('fullpath');
    repo_dir = fileparts(fileparts(fileparts(this_file)));
end
if nargin < 2 || isempty(molecule)
    molecule = 'sucrose';
end

config_file = fullfile(repo_dir, 'data', molecule, [molecule '_config.json']);
if ~isfile(config_file)
    error('Missing carbohydrate configuration: %s', config_file);
end
config = jsondecode(fileread(config_file));
if ~isfield(config, 'name') || ~strcmp(config.name, molecule)
    error('Configuration %s does not declare the expected molecule %s.', config_file, molecule);
end
end
