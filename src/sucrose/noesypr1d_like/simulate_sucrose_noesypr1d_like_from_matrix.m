% simulate_sucrose_noesypr1d_like_from_matrix.m
%
% Reads a GISSMO/A2/A5-style spin matrix file:
%   diagonal     = 1H chemical shifts in ppm
%   off-diagonal = J couplings in Hz
%
% Then simulates a Bruker noesypr1d-like 1D 1H spectrum using
% the real Bruker acquisition parameters from the 1100 MHz DATA_NOESY file.

close all;
clc;

script_dir = fileparts(mfilename('fullpath'));
repo_root = getenv('EDISON_REPO_ROOT');
if isempty(repo_root)
    repo_root = fileparts(fileparts(fileparts(script_dir)));
end
addpath(fullfile(repo_root, 'src', 'common'));
config = load_carbohydrate_config(repo_root, 'sucrose');
legacy = config.legacy_noesypr1d_like;
sequence = config.sequence;

% ------------------------------------------------------------
% User-selectable matrix file
% ------------------------------------------------------------
if ~exist('matrix_file', 'var')
    matrix_file = getenv('SUCROSE_MATRIX_FILE');
end

if isempty(matrix_file)
    script_dir_default = fileparts(mfilename('fullpath'));
    repo_dir_default = fileparts(fileparts(fileparts(script_dir_default)));
    matrix_file = fullfile(repo_dir_default, 'data', 'sucrose', 'matrix', ...
        'sucrose_spin_matrix_GISSMO_14x14.txt');
end

fprintf('Reading spin matrix from: %s\n', matrix_file);

M = readmatrix(matrix_file, 'FileType', 'text');

% Remove fully empty rows/columns if readmatrix picked up blanks
M = M(any(~isnan(M), 2), :);
M = M(:, any(~isnan(M), 1));

if size(M,1) ~= size(M,2)
    error('Spin matrix must be square. Read size was %d x %d.', size(M,1), size(M,2));
end

nspins = size(M,1);
fprintf('Matrix size: %d x %d\n', nspins, nspins);

shifts = diag(M);

if any(isnan(shifts))
    error('Matrix diagonal contains NaN values. Diagonal must contain chemical shifts in ppm.');
end

% Build a clean J matrix without double-counting.
% If both triangles contain the same J, use it once.
% If only one triangle contains a J, use the nonzero value.
J = zeros(nspins);

for a = 1:nspins
    for b = a+1:nspins
        jab = M(a,b);
        jba = M(b,a);

        if isnan(jab); jab = 0; end
        if isnan(jba); jba = 0; end

        if abs(jab) > 1e-12 && abs(jba) > 1e-12
            J(a,b) = 0.5 * (jab + jba);
        elseif abs(jab) > 1e-12
            J(a,b) = jab;
        elseif abs(jba) > 1e-12
            J(a,b) = jba;
        end
    end
end

fprintf('Chemical shift range: %.4f to %.4f ppm\n', min(shifts), max(shifts));
fprintf('Number of nonzero J couplings used: %d\n', nnz(triu(J,1)));

% ------------------------------------------------------------
% Output names
% ------------------------------------------------------------
[~, matrix_stem, ~] = fileparts(matrix_file);
matrix_stem = regexprep(matrix_stem, '[^A-Za-z0-9_]+', '_');

% ------------------------------------------------------------
% Line broadening, adjustable from shell:
%   LB_HZ=0.5 matlab -batch "run(...)"
% ------------------------------------------------------------
if ~exist('lb', 'var')
    lb_env = getenv('LB_HZ');
    if isempty(lb_env)
        lb = legacy.default_lb_hz;
    else
        lb = str2double(lb_env);
    end
end

if isnan(lb) || lb <= 0
    error('LB_HZ must be a positive number.');
end

lb_tag = sprintf('lb%.2f', lb);
lb_tag = strrep(lb_tag, '.', 'p');

fprintf('Line broadening LB = %.3f Hz\n', lb);

out_dir = fullfile(repo_root, 'outputs', 'sucrose', 'noesypr1d_like_from_matrix');

if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

out_prefix = fullfile(out_dir, ['noesypr1d_like_' matrix_stem '_' lb_tag]);

fprintf('Output prefix: %s\n', out_prefix);

% ------------------------------------------------------------
% Add Spinach
% ------------------------------------------------------------
spinach_root = getenv('SPINACH_ROOT');
if isempty(spinach_root)
    error('Set SPINACH_ROOT before running the legacy noesypr1d-like script.');
end
addpath(genpath(spinach_root));

% ------------------------------------------------------------
% Bruker acquisition parameters from environment variables
% ------------------------------------------------------------
% Required: SFO1_MHZ, O1_HZ, SW_HZ, TD_POINTS. Optional: BF1_MHZ.
SFO1_env = getenv('SFO1_MHZ');
O1_env = getenv('O1_HZ');
SW_env = getenv('SW_HZ');
TD_env = getenv('TD_POINTS');
BF1_env = getenv('BF1_MHZ');
if any(cellfun(@isempty, {SFO1_env, O1_env, SW_env, TD_env}))
    error('Need SFO1_MHZ, O1_HZ, SW_HZ, and TD_POINTS environment variables.');
end
SFO1_MHz = str2double(SFO1_env);
O1_Hz = str2double(O1_env);
SW_h_Hz = str2double(SW_env);
TD_real = str2double(TD_env);
if isempty(BF1_env), BF1_MHz = SFO1_MHz; else, BF1_MHz = str2double(BF1_env); end
if any(isnan([SFO1_MHz O1_Hz SW_h_Hz TD_real BF1_MHz]))
    error('One or more acquisition environment variables is not numeric.');
end
n_complex = TD_real / 2;
SW_ppm = SW_h_Hz / SFO1_MHz;

center_ppm = O1_Hz / BF1_MHz;
width_ppm = SW_ppm;

fprintf('Bruker center ppm from O1/BF1 = %.6f ppm\n', center_ppm);
fprintf('Bruker spectral width = %.6f ppm\n', width_ppm);

% ------------------------------------------------------------
% Spin system
% ------------------------------------------------------------
sys.magnet = SFO1_MHz / config.spinach.proton_gamma_mhz_per_t;
sys.isotopes = repmat({'1H'}, 1, nspins);

labels = cell(1, nspins);
for k = 1:nspins
    labels{k} = sprintf('Spin_%02d', k);
end
sys.labels = labels;

% Chemical shifts
inter.zeeman.scalar = cell(1, nspins);
for k = 1:nspins
    inter.zeeman.scalar{k} = shifts(k);
end

% J couplings
inter.coupling.scalar = cell(nspins);
for a = 1:nspins
    for b = a+1:nspins
        if abs(J(a,b)) > 1e-12
            inter.coupling.scalar{a,b} = J(a,b);
        end
    end
end

% ------------------------------------------------------------
% Approximate relaxation model
% ------------------------------------------------------------
use_relaxation = true;

if use_relaxation
    inter.relaxation = {'t1_t2'};
    inter.rlx_keep = 'secular';
    inter.equilibrium = 'zero';

    % Cell arrays are required by this Spinach version.
    inter.r1_rates = num2cell(legacy.r1_hz * ones(1, nspins));
    inter.r2_rates = num2cell(legacy.r2_hz * ones(1, nspins));
end

% ------------------------------------------------------------
% Basis set
% ------------------------------------------------------------
bas.formalism = 'sphten-liouv';
bas.approximation = config.spinach.basis_approximation;
bas.connectivity = config.spinach.basis_connectivity;
bas.space_level = config.spinach.basis_space_level;

% ------------------------------------------------------------
% Build spin system
% ------------------------------------------------------------
spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

% ------------------------------------------------------------
% NOESYPR1D-like sequence parameters
% ------------------------------------------------------------
parameters.spins = {'1H'};

parameters.offset = O1_Hz;
parameters.sweep = SW_h_Hz;
parameters.npoints = n_complex;
parameters.zerofill = legacy.zerofill_points;
parameters.axis_units = 'ppm';
parameters.invert_axis = 1;

% Real Bruker delays from acqus:
% D1 = 2 s
% D8 = 0.05 s
parameters.d1 = sequence.d1_s;
parameters.tmix = sequence.tmix_s;

parameters.rho0 = state(spin_system, 'Lz', '1H');
parameters.coil = state(spin_system, 'L+', '1H');

parameters.decouple = {};
parameters.verbose = 1;

% ------------------------------------------------------------
% Simulate FID
% ------------------------------------------------------------
fid = liquid(spin_system, @noesypr1d_like, parameters, 'nmr');
fid = fid(:);

% ------------------------------------------------------------
% Process FID
% ------------------------------------------------------------
% lb is defined near the top of the script from LB_HZ or default 1.0 Hz
dt = 1 / parameters.sweep;
t = (0:numel(fid)-1).' * dt;

fid_apod = fid .* exp(-pi * lb * t);
spec = real(fftshift(fft(fid_apod, parameters.zerofill)));

% Auto-flip if mostly negative
if abs(min(spec)) > abs(max(spec))
    spec = -spec;
end

ppm_axis = linspace(center_ppm - width_ppm/2, ...
                   center_ppm + width_ppm/2, ...
                   parameters.zerofill).';

mx = max(abs(spec));
if mx == 0 || isnan(mx)
    error('Simulated spectrum is zero or NaN before normalization.');
end
spec = spec / mx;

% Save full simulated spectrum
full_file = [out_prefix '_full_simulated_spectrum.txt'];
writematrix([ppm_axis, spec], full_file, 'Delimiter', 'tab');
fprintf('Saved %s\n', full_file);

% Save sucrose-region crop
crop_min = legacy.crop_region_ppm(1);
crop_max = legacy.crop_region_ppm(2);
crop_mask = ppm_axis >= crop_min & ppm_axis <= crop_max;

crop_file = [out_prefix '_sucrose_region_simulated_spectrum.txt'];
writematrix([ppm_axis(crop_mask), spec(crop_mask)], crop_file, 'Delimiter', 'tab');
fprintf('Saved %s\n', crop_file);

% ------------------------------------------------------------
% Compare to experimental DATA_NOESY
% ------------------------------------------------------------
exp_file = getenv('EXP_FILE');

if isfile(exp_file)
    E = readmatrix(exp_file);
    E = E(:, 1:min(2, size(E,2)));
    E = E(all(~isnan(E), 2), :);

    exp_ppm = E(:,1);
    exp_y = E(:,2);

    exp_y = exp_y - median(exp_y);
    exp_y = exp_y / max(abs(exp_y));

    [sim_ppm_sorted, sim_idx] = sort(ppm_axis);
    sim_y_sorted = spec(sim_idx);

    [exp_ppm_sorted, exp_idx] = sort(exp_ppm);
    exp_y_sorted = exp_y(exp_idx);

    sim_on_exp = interp1(sim_ppm_sorted, sim_y_sorted, exp_ppm_sorted, 'pchip', 0);

    comparison = [exp_ppm_sorted, exp_y_sorted, sim_on_exp, exp_y_sorted - sim_on_exp];

    comparison_file = [out_prefix '_vs_experiment_comparison.txt'];
    writematrix(comparison, comparison_file, 'Delimiter', 'tab');
    fprintf('Saved %s\n', comparison_file);

    % Shape metrics over sucrose region
    region_mask = exp_ppm_sorted >= legacy.comparison_region_ppm(1) & ...
        exp_ppm_sorted <= legacy.comparison_region_ppm(2);
    e = exp_y_sorted(region_mask);
    s = sim_on_exp(region_mask);

    rmse = sqrt(mean((e - s).^2));
    mae = mean(abs(e - s));
    cos_sim = dot(e, s) / (norm(e) * norm(s));
    cos_loss = 1 - cos_sim;

    metrics_file = [out_prefix '_shape_metrics.txt'];
    fidm = fopen(metrics_file, 'w');
    fprintf(fidm, 'matrix_file\t%s\n', matrix_file);
    fprintf(fidm, 'nspins\t%d\n', nspins);
    fprintf(fidm, 'D1_seconds\t%.6f\n', parameters.d1);
    fprintf(fidm, 'D8_tmix_seconds\t%.6f\n', parameters.tmix);
    fprintf(fidm, 'rmse_3p35_to_5p45\t%.10f\n', rmse);
    fprintf(fidm, 'mae_3p35_to_5p45\t%.10f\n', mae);
    fprintf(fidm, 'cosine_similarity_3p35_to_5p45\t%.10f\n', cos_sim);
    fprintf(fidm, 'cosine_loss_3p35_to_5p45\t%.10f\n', cos_loss);
    fclose(fidm);
    fprintf('Saved %s\n', metrics_file);

    % Plot full sucrose region
    fig = figure('Visible', 'off');
    plot(exp_ppm_sorted, exp_y_sorted, 'LineWidth', 1.0);
    hold on;
    plot(exp_ppm_sorted, sim_on_exp, 'LineWidth', 1.0);
    set(gca, 'XDir', 'reverse');
    xlabel('1H chemical shift / ppm');
    ylabel('normalized intensity');
    title(['DATA NOESY vs noesypr1d-like simulation: ' matrix_stem], 'Interpreter', 'none');
    legend('experiment', 'simulation');
    grid on;

    plot_file = [out_prefix '_vs_experiment.png'];
    saveas(fig, plot_file);
    close(fig);
    fprintf('Saved %s\n', plot_file);

    % Plot lower sugar region
    lower_mask = exp_ppm_sorted >= legacy.lower_region_ppm(1) & ...
        exp_ppm_sorted <= legacy.lower_region_ppm(2);
    fig = figure('Visible', 'off');
    plot(exp_ppm_sorted(lower_mask), exp_y_sorted(lower_mask), 'LineWidth', 1.0);
    hold on;
    plot(exp_ppm_sorted(lower_mask), sim_on_exp(lower_mask), 'LineWidth', 1.0);
    set(gca, 'XDir', 'reverse');
    xlabel('1H chemical shift / ppm');
    ylabel('normalized intensity');
    title(['Lower sugar region: ' matrix_stem], 'Interpreter', 'none');
    legend('experiment', 'simulation');
    grid on;

    lower_plot_file = [out_prefix '_lower_3p80_3p35.png'];
    saveas(fig, lower_plot_file);
    close(fig);
    fprintf('Saved %s\n', lower_plot_file);

else
    fprintf('Experimental file not found: %s\n', exp_file);
end

fprintf('Simulation complete.\n');

% ========================================================================
% Local custom pulse sequence: 1D NOESYPR1D-like
% ========================================================================
function fid = noesypr1d_like(spin_system, parameters, H, R, K)

    L = H + 1i*R + 1i*K;

    coil = parameters.coil;

    Lp = operator(spin_system, 'L+', parameters.spins{1});
    Lx = (Lp + Lp') / 2;
    Ly = (Lp - Lp') / (2i);

    rho = parameters.rho0;

    % Relaxation delay / approximate presaturation period.
    if isfield(parameters, 'd1') && parameters.d1 > 0
        rho = evolution(spin_system, 1i*R + 1i*K, [], ...
                        rho, parameters.d1, 1, 'final');
    end

    % Approximate 1D NOESY-like preparation:
    % 90x - 90x - spoil - mix - spoil - 90y - acquire
    rho = step(spin_system, Lx, rho, pi/2);
    rho = step(spin_system, Lx, rho, pi/2);

    rho = homospoil(spin_system, rho, 'destroy');

    if isfield(parameters, 'tmix') && parameters.tmix > 0
        rho = evolution(spin_system, 1i*R + 1i*K, [], ...
                        rho, parameters.tmix, 1, 'final');
    end

    rho = homospoil(spin_system, rho, 'destroy');

    rho = step(spin_system, Ly, rho, pi/2);

    timestep = 1 / parameters.sweep;

    fid = evolution(spin_system, L, coil, ...
                    rho, timestep, parameters.npoints - 1, ...
                    'observable');
end
