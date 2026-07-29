function results = run_carbohydrate_spinach_field(molecule, field_label, matrix_override)
% Run the shared Spinach forward model for one molecule and Bruker field.
% The molecule-specific matrix, blocks, and acquisition metadata come from
% data/<molecule> and outputs/<molecule>/prepared.

if nargin < 1 || isempty(molecule), molecule = 'sucrose'; end
if nargin < 2 || isempty(field_label), error('Provide a field label, e.g. ''700''.'); end
if nargin < 3, matrix_override = ''; end

common_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(common_dir));
addpath(common_dir);
addpath(fullfile(repo_dir, 'src', 'sucrose'));

config = load_carbohydrate_config(repo_dir, molecule);
metadata_file = fullfile(repo_dir, 'outputs', molecule, 'prepared', ...
    'preparation_summary.csv');
if ~isfile(metadata_file)
    error('Missing %s. Run prepare_carbohydrate_spectra.py first.', metadata_file);
end
metadata = readtable(metadata_file, 'TextType', 'string');
field_mhz = str2double(string(field_label));
row = find(metadata.field_mhz == field_mhz, 1, 'first');
if isempty(row)
    error('No metadata row for %s MHz in %s.', field_label, metadata_file);
end

spinach_root = getenv('SPINACH_ROOT');
if isempty(spinach_root)
    error('Set SPINACH_ROOT before running Spinach.');
end

p.molecule = molecule;
if isempty(matrix_override)
    matrix_tag = '';
else
    matrix_tag = '_candidate';
end
p.project_dir = fullfile(repo_dir, 'outputs', molecule, ...
    [field_label 'MHz_spinach' matrix_tag]);
p.spinach_root = spinach_root;
if ~exist(p.project_dir, 'dir'), mkdir(p.project_dir); end

p.matrix_file = char(string(config.matrix_file));
if ~isempty(matrix_override)
    p.matrix_file = char(string(matrix_override));
end
% Always pass an absolute matrix path. Spinach may change MATLAB's current
% folder while creating its system, so a relative path can work during setup
% and then fail inside simulate_sucrose_spinach_fft.
is_absolute_matrix_path = startsWith(p.matrix_file, '/') || ...
    startsWith(p.matrix_file, char(92)) || ...
    (numel(p.matrix_file) >= 2 && p.matrix_file(2) == ':');
if ~is_absolute_matrix_path
    p.matrix_file = fullfile(repo_dir, p.matrix_file);
end
if ~isfile(p.matrix_file)
    error('Configured matrix file does not exist: %s', p.matrix_file);
end

% If the BMRB/GISSMO seed includes deposited field-specific spectra, pass the
% matching JSON into the Spinach diagnostic so the overlay contains all three
% curves: experiment, current Spinach run, and published GISSMO simulation.
if isfield(config, 'gissmo_simulation_dir')
    gissmo_dir = char(string(config.gissmo_simulation_dir));
    if ~(startsWith(gissmo_dir, '/') || startsWith(gissmo_dir, char(92)) || ...
            (numel(gissmo_dir) >= 2 && gissmo_dir(2) == ':'))
        gissmo_dir = fullfile(repo_dir, gissmo_dir);
    end
    gissmo_file = fullfile(gissmo_dir, sprintf('sim_%dMHz.json', round(field_mhz)));
    if isfile(gissmo_file)
        p.gissmo_sim_file = gissmo_file;
    else
        warning('No deposited GISSMO spectrum found for %g MHz: %s', field_mhz, gissmo_file);
    end
end
p.atom_ids = config.atom_ids;
% jsondecode returns a numeric matrix for equal-length blocks and a cell
% array for ragged blocks (for example sucrose's 7-, 5-, and 2-spin
% components). Normalize both forms to a flat cell array of numeric indices.
if iscell(config.blocks)
    p.blocks = config.blocks;
else
    p.blocks = num2cell(config.blocks, 2);
end
p.sample_label = molecule;
p.field_label = [field_label ' MHz'];
p.plot_prefix = [molecule '_' field_label 'MHz_spinach' matrix_tag];
p.suppress_sucrose_diagnostics = true;

p.SFO1_MHz = metadata.sfo1_mhz(row);
p.O1_Hz = metadata.o1_hz(row);
p.SW_Hz = metadata.sw_acq_hz(row);
p.ncomplex = metadata.ncomplex(row);
p.TD = metadata.points(row);
p.exp_file = fullfile(repo_dir, 'data', molecule, ...
    char(metadata.relative_dir(row)), 'pdata', char(string(metadata.procno(row))), '1r');
p.EXP_SI = metadata.points(row);
p.EXP_SF = metadata.sf_mhz(row);
p.EXP_SW_p = metadata.sw_hz(row);
p.EXP_OFFSET = metadata.offset_dss_ppm(row);
p.EXP_BYTORDP = metadata.bytordp(row);
p.EXP_DTYPP = metadata.dtypp(row);
p.EXP_NC_proc = metadata.nc_proc(row);

p.water_win = config.processing.water_region_ppm;
p.artifact_win = config.processing.artifact_region_ppm;
p.anomeric_win = config.processing.anomeric_region_ppm;
p.crowded_win = config.processing.crowded_region_ppm;
p.sucrose_region = config.processing.fit_region_ppm;
p.bas_approximation = 'IK-2';
p.bas_connectivity = 'scalar_couplings';
p.bas_space_level = 3;
p.parallel_workers = 1;
p.lb_Hz_guess = 1.0;

p.acquisition_mode = 'direct';
p.noesy_sequence_model = 'homospoil';
p.fit_receiver_phase = false;
p.noesy_tmix_s = 0.05;
p.noesy_d1_s = 0;
p.lineshape_model = 'voigt';
if ismember('pulse_program', metadata.Properties.VariableNames)
    pulse = lower(char(metadata.pulse_program(row)));
    if contains(pulse, 'noesy')
        p.acquisition_mode = 'noesypr1d';
        p.fit_receiver_phase = true;
        p.noesy_d1_s = 2.0;
    end
end

results = simulate_sucrose_spinach_fft(p);

% Add runner-level provenance fields so the all-field aggregator can collect
% results without reopening the JSON summary.
results.field_mhz = field_mhz;
results.field_label = p.field_label;
results.acquisition_mode = p.acquisition_mode;
results.pulse_program = char(metadata.pulse_program(row));
results.ppm_offset_fitted = results.exp_ppm_offset_fitted;

% Persist only scalar provenance/fit results. The large FID and spectral
% arrays remain in MATLAB's return value and are not written into JSON.
summary = struct();
summary.molecule = molecule;
summary.field_mhz = field_mhz;
summary.field_label = p.field_label;
summary.matrix_file = p.matrix_file;
summary.nspins = results.nspins;
summary.blocks = p.blocks;
summary.acquisition_mode = p.acquisition_mode;
summary.pulse_program = char(metadata.pulse_program(row));
summary.r_spinach_vs_expt = results.r_spinach_vs_expt;
summary.rmse_spinach_vs_expt = results.rmse_spinach_vs_expt;
summary.lbL_Hz = results.lbL_Hz;
summary.lbG_Hz = results.lbG_Hz;
summary.receiver_phase_deg = results.receiver_phase_deg;
summary.ppm_offset_fitted = results.exp_ppm_offset_fitted;
summary_file = fullfile(p.project_dir, ...
    [molecule '_' field_label 'MHz_spinach_summary.json']);
fid = fopen(summary_file, 'w');
if fid < 0, error('Unable to write Spinach summary: %s', summary_file); end
fwrite(fid, jsonencode(summary), 'char');
fwrite(fid, newline, 'char');
fclose(fid);
results.summary_file = summary_file;
fprintf('Wrote %s\n', summary_file);
end
