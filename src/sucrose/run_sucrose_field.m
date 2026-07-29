function results = run_sucrose_field(field_label, make_plots)
% Configure one official sucrose field and call the shared Spinach engine.
% Primary set: 600/900/1100 MHz. 800 MHz is available as a separate
% condition and can be run for robustness without changing the matrix.

if nargin < 2 || isempty(make_plots)
    make_plots = true;
end

src_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(src_dir));
addpath(src_dir);
addpath(fullfile(repo_dir, 'src', 'common'));

config = load_carbohydrate_config(repo_dir, 'sucrose');
sequence = config.sequence;
spinach = config.spinach;

spinach_root = getenv('SPINACH_ROOT');
% Prefer an explicitly configured checkout, but make the repository
% self-contained for students by falling back to the bundled library.
if isempty(spinach_root) || exist(spinach_root, 'dir') ~= 7
    bundled_spinach = fullfile(repo_dir, 'lib', 'Spinach-2.10.1');
    if exist(bundled_spinach, 'dir') == 7
        spinach_root = bundled_spinach;
        fprintf('SPINACH_ROOT not set; using bundled Spinach: %s\n', spinach_root);
    else
        error(['Spinach not found. Set SPINACH_ROOT or install it at ' ...
            '%s.'], bundled_spinach);
    end
end

p.project_dir = fullfile(repo_dir, 'outputs', 'sucrose', [field_label 'MHz']);
p.spinach_root = spinach_root;
if ~exist(p.project_dir, 'dir'), mkdir(p.project_dir); end

seed_manifest_file = fullfile(repo_dir, 'outputs', 'sucrose', 'seed_selection.json');
if isfile(seed_manifest_file)
    seed_manifest = jsondecode(fileread(seed_manifest_file));
    allow_provisional = isfield(config, 'seed_selection') && ...
        isfield(config.seed_selection, 'allow_provisional') && config.seed_selection.allow_provisional;
    if ~strcmp(seed_manifest.status, 'READY') && ...
            ~(allow_provisional && strcmp(seed_manifest.status, 'PROVISIONAL_SEED'))
        error('Seed selection is %s; provide a numeric matrix and spectral provenance before fitting.', seed_manifest.status);
    elseif strcmp(seed_manifest.status, 'PROVISIONAL_SEED')
        fprintf('WARNING: fitting a provisional Bubb+spectra seed; this is not yet a publishable GISSMO matrix.\n');
    end
    p.matrix_file = char(seed_manifest.matrix_file);
    % Seed manifests intentionally store repository-relative paths.  Resolve
    % them before simulate_sucrose_spinach_fft changes MATLAB's current
    % folder to the field-specific output directory.
    if ~startsWith(p.matrix_file, filesep)
        p.matrix_file = fullfile(repo_dir, p.matrix_file);
    end
else
p.matrix_file = fullfile(repo_dir, config.matrix_file);
end
p.atom_ids = config.atom_ids;
% MATLAB jsondecode returns a cell array when the JSON block rows have
% different lengths (as in sucrose's 7/5/2-spin decomposition).  Preserve
% those numeric row vectors; wrapping an already-cell array with num2cell
% creates nested cells and later makes idx(k) invalid in the Spinach engine.
if iscell(config.blocks)
    p.blocks = config.blocks;
else
    p.blocks = num2cell(config.blocks, 2);
end
p.acquisition_mode = char(sequence.acquisition_mode);
p.noesy_sequence_model = char(sequence.sequence_model);
p.noesy_tmix_s = sequence.tmix_s;
p.noesy_d1_s = sequence.d1_s;
p.fit_receiver_phase = sequence.fit_receiver_phase;
p.lineshape_model = char(sequence.lineshape_model);
p.lb_Hz_guess = sequence.lb_hz_guess;
p.water_win = config.processing.water_region_ppm;
p.artifact_win = config.processing.artifact_region_ppm;
p.anomeric_win = config.processing.anomeric_region_ppm;
p.crowded_win = config.processing.crowded_region_ppm;
p.bas_approximation = char(spinach.basis_approximation);
p.bas_connectivity = char(spinach.basis_connectivity);
p.bas_space_level = spinach.basis_space_level;
p.parallel_workers = sequence.parallel_workers;

% Acquisition and processing numbers are generated from the Bruker acqus,
% procs, and 1r files by prepare_sucrose_spectra.py.  This keeps field
% metadata out of the driver and makes the same pattern reusable for other
% carbohydrates.
metadata_file = fullfile(repo_dir, 'outputs', 'sucrose', 'prepared', ...
    'preparation_summary.csv');
if ~isfile(metadata_file)
    error('Missing %s. Run prepare_sucrose_spectra.py first.', metadata_file);
end
metadata = readtable(metadata_file, 'TextType', 'string');
field_mhz = str2double(field_label);
row = find(metadata.field_mhz == field_mhz, 1, 'first');
if isempty(row)
    error('No metadata row for %s MHz in %s.', field_label, metadata_file);
end

p.SFO1_MHz = metadata.sfo1_mhz(row);
p.O1_Hz = metadata.o1_hz(row);
p.SW_Hz = metadata.sw_acq_hz(row);
p.ncomplex = metadata.ncomplex(row);
p.TD = metadata.points(row);
p.exp_file = fullfile(repo_dir, 'data', 'sucrose', ...
    char(metadata.relative_dir(row)), 'pdata', ...
    char(string(metadata.procno(row))), '1r');
p.EXP_SI = metadata.points(row);
p.EXP_SF = metadata.sf_mhz(row);
p.EXP_SW_p = metadata.sw_hz(row);
p.EXP_OFFSET = metadata.offset_dss_ppm(row);
p.EXP_BYTORDP = metadata.bytordp(row);
p.EXP_DTYPP = metadata.dtypp(row);
p.EXP_NC_proc = metadata.nc_proc(row);

p.field_label = [field_label ' MHz'];
p.plot_prefix = ['sucrose_' field_label 'MHz'];
p.make_plots = logical(make_plots);
results = simulate_sucrose_spinach_fft(p);
end
