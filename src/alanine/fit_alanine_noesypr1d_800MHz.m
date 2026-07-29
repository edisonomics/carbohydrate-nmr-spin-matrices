% fit_alanine_noesypr1d_experiment.m
%
% Fit the real 800 MHz alanine noesypr1d processed spectrum using the
% Bruker-style noesypr1d pulse block validated in diagnose_alanine_noesypr1d.

clear; close all; clc;

src_dir      = fileparts(mfilename('fullpath'));
repo_dir     = fileparts(fileparts(src_dir));
project_dir  = fullfile(repo_dir, 'outputs', 'alanine', '800MHz');
spinach_root = getenv('SPINACH_ROOT');
exp_1r       = fullfile(repo_dir, 'data', 'alanine', '800_MHz', '11', 'pdata', '1', '1r');

if isempty(spinach_root), error('Set SPINACH_ROOT to the Spinach installation.'); end
if ~exist(project_dir, 'dir'), mkdir(project_dir); end
addpath(src_dir);

cd(project_dir);

if exist(spinach_root, 'dir') && exist('create', 'file') ~= 2
    addpath(genpath(spinach_root));
end
if exist('create', 'file') ~= 2
    error('Spinach create() not found. Check spinach_root.');
end
if ~isfile(exp_1r)
    error('Experimental 1r not found: %s', exp_1r);
end

%% Acquisition settings from data/alanine/800_MHz/11
SFO1_MHz = 799.713758637;
O1_Hz    = 3758.637;
SW_Hz    = 9615.385;
TD       = 32768;
ncomplex = TD/2;
tmix_s   = 0.05;

methyl_win = [1.30 1.65];
alpha_win  = [3.55 3.95];

fprintf('\n===== Alanine 800 MHz noesypr1d experiment fit =====\n');
fprintf('Experimental file: %s\n', exp_1r);
fprintf('SFO1 %.9f MHz, O1 %.3f Hz, SW %.1f Hz, TD %d, d8 %.3f s\n', ...
    SFO1_MHz, O1_Hz, SW_Hz, TD, tmix_s);

%% Build alanine AX3 spin system
B0_T = SFO1_MHz / 42.57747892;

sys.magnet   = B0_T;
sys.isotopes = {'1H', '1H', '1H', '1H'};

delta_Ha_ppm  = 3.7680;
delta_CH3_ppm = 1.4655;
J_Ha_CH3_Hz   = 7.234;

inter.zeeman.scalar = {delta_Ha_ppm, delta_CH3_ppm, delta_CH3_ppm, delta_CH3_ppm};
inter.coupling.scalar      = cell(4, 4);
inter.coupling.scalar{1,2} = J_Ha_CH3_Hz;
inter.coupling.scalar{1,3} = J_Ha_CH3_Hz;
inter.coupling.scalar{1,4} = J_Ha_CH3_Hz;

bas.formalism     = 'sphten-liouv';
bas.approximation = 'none';

spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

parameters.spins       = {'1H'};
parameters.offset      = O1_Hz;
parameters.sweep       = SW_Hz;
parameters.npoints     = ncomplex;
parameters.zerofill    = TD;
parameters.axis_units  = 'ppm';
parameters.invert_axis = 1;
parameters.decouple    = {};
parameters.rho0        = state(spin_system, 'Lz', '1H', 'cheap');
parameters.coil        = state(spin_system, 'L+', '1H', 'cheap');
parameters.tmix        = tmix_s;
parameters.pulse_sign  = 1;
parameters.receiver_sign = 1;
parameters.phase_cycle = true;
parameters.crusher_mode = 'none';   % real noesypr1d has no gradients

fprintf('\nRunning Spinach noesypr1d pulse block once...\n');
tic;
fid_noesy = liquid(spin_system, @alanine_noesypr1d_acquire, parameters, 'nmr');
fid_noesy = fid_noesy(:);
fprintf('NOESY FID done in %.2f s; max |FID| %.6g\n', toc, max(abs(fid_noesy)));

%% Read and normalise experiment
[exp_ppm, exp_y] = read_bruker_1r(exp_1r);
exp_norm = normalize_trace(exp_y, in_window(exp_ppm, methyl_win));

%% Native ppm axis and fitting masks
[ppm_axis, ~] = process_fid(fid_noesy, SFO1_MHz, O1_Hz, SW_Hz, TD, 1.5);
score_mask = in_window(ppm_axis, methyl_win) | in_window(ppm_axis, alpha_win);
center_ppm = mean(ppm_axis(score_mask));

%% Fit LB, ppm offset, and zero-order receiver phase.
% Scale and a linear baseline are solved analytically inside the objective.
q0 = [log(1.4), -0.002, pi/2];
opts = optimset('Display', 'iter', 'TolX', 1e-7, 'TolFun', 1e-8, ...
    'MaxFunEvals', 600, 'MaxIter', 300);
objective = @(q) fit_rmse(q, fid_noesy, SFO1_MHz, O1_Hz, SW_Hz, TD, ...
    ppm_axis, exp_ppm, exp_norm, score_mask, center_ppm);

fprintf('\nFitting LB, ppm offset, and receiver phase...\n');
q_opt = fminsearch(objective, q0, opts);

[rmse_fit, fit_details] = fit_rmse(q_opt, fid_noesy, SFO1_MHz, O1_Hz, SW_Hz, TD, ...
    ppm_axis, exp_ppm, exp_norm, score_mask, center_ppm);

lb_Hz      = exp(q_opt(1));
ppm_offset = q_opt(2);
phase_deg  = wrap_to_180(q_opt(3) * 180/pi);

r_fit = trace_corr(fit_details.target, fit_details.fit);

fprintf('\n===== Fitted noesypr1d -> experiment =====\n');
fprintf('LB          : %.6f Hz\n', lb_Hz);
fprintf('ppm offset  : %+ .7f ppm  (added to experimental axis)\n', ppm_offset);
fprintf('phase       : %+ .3f deg\n', phase_deg);
fprintf('scale       : %.8g\n', fit_details.coef(1));
fprintf('baseline b0 : %.8g\n', fit_details.coef(2));
fprintf('baseline b1 : %.8g per ppm\n', fit_details.coef(3));
fprintf('r           : %.6f\n', r_fit);
fprintf('RMSE        : %.8f\n', rmse_fit);

%% Build final full-grid traces
[~, model_complex] = process_fid(fid_noesy, SFO1_MHz, O1_Hz, SW_Hz, TD, lb_Hz);
model_real = real(exp(1i*q_opt(3)) * model_complex);
exp_on_model = interp1(exp_ppm + ppm_offset, exp_norm, ppm_axis, 'linear', 0);

X_full = [model_real(:), ones(numel(ppm_axis), 1), ppm_axis(:) - center_ppm];
model_fit = X_full * fit_details.coef;
residual = model_fit - exp_on_model;

%% Save numeric overlay
overlay = table(ppm_axis(:), exp_on_model(:), model_fit(:), residual(:), ...
    'VariableNames', {'ppm', 'experiment_norm_interp', 'noesypr1d_fit', 'residual_model_minus_expt'});
writetable(overlay, fullfile(project_dir, 'alanine_noesypr1d_experiment_fit_overlay.csv'));

summary = table(lb_Hz, ppm_offset, phase_deg, fit_details.coef(1), ...
    fit_details.coef(2), fit_details.coef(3), r_fit, rmse_fit, ...
    'VariableNames', {'lb_Hz', 'ppm_offset', 'phase_deg', 'scale', ...
    'baseline_b0', 'baseline_b1', 'r_noesypr1d_vs_expt', 'rmse_noesypr1d_vs_expt'});
writetable(summary, fullfile(project_dir, 'alanine_noesypr1d_experiment_fit_summary.csv'));

save(fullfile(project_dir, 'alanine_noesypr1d_experiment_fit.mat'), ...
    'summary', 'overlay', 'fid_noesy', 'ppm_axis', 'exp_ppm', 'exp_norm', ...
    'model_fit', 'exp_on_model', 'residual', 'fit_details', 'q_opt');

%% Plot
fig = figure('Color', 'w', 'InvertHardcopy', 'off', 'Position', [70 70 1700 960]);
tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Alanine 800 MHz noesypr1d fit to experiment | r %.4f, RMSE %.5f', ...
    r_fit, rmse_fit), 'Color', 'k', 'FontSize', 19, 'FontWeight', 'bold');

nexttile([1 2]);
hold on;
h1 = plot(ppm_axis, exp_on_model, '-', 'Color', [0 0 0], 'LineWidth', 2.2);
h2 = plot(ppm_axis, model_fit, '--', 'Color', [0.86 0.05 0.05], 'LineWidth', 2.6);
style_axis(gca, 18);
xlim([1.1 4.1]);
ylim([-0.15 1.15]);
xlabel('^{1}H chemical shift (ppm)', 'Color', 'k', 'FontWeight', 'bold');
ylabel('normalised intensity', 'Color', 'k', 'FontWeight', 'bold');
title('Full alanine window', 'Color', 'k', 'FontWeight', 'bold');
subtitle(sprintf('LB %.3f Hz, offset %+ .5f ppm, phase %+ .1f deg', ...
    lb_Hz, ppm_offset, phase_deg), 'Color', 'k', 'FontSize', 13, 'FontWeight', 'bold');
lg = legend([h1 h2], 'experiment', 'Spinach noesypr1d fit', 'Location', 'northwest');
style_legend(lg);

nexttile;
hold on;
plot(ppm_axis, exp_on_model, '-', 'Color', [0 0 0], 'LineWidth', 2.3);
plot(ppm_axis, model_fit, '--', 'Color', [0.86 0.05 0.05], 'LineWidth', 2.6);
style_axis(gca, 18);
xlim(methyl_win);
ylim([-0.10 1.15]);
xlabel('ppm', 'Color', 'k', 'FontWeight', 'bold');
ylabel('normalised intensity', 'Color', 'k', 'FontWeight', 'bold');
title('Methyl zoom', 'Color', 'k', 'FontWeight', 'bold');

nexttile;
hold on;
plot(ppm_axis, exp_on_model, '-', 'Color', [0 0 0], 'LineWidth', 2.3);
plot(ppm_axis, model_fit, '--', 'Color', [0.86 0.05 0.05], 'LineWidth', 2.6);
style_axis(gca, 18);
xlim(alpha_win);
ylim([-0.05 0.35]);
xlabel('ppm', 'Color', 'k', 'FontWeight', 'bold');
ylabel('normalised intensity', 'Color', 'k', 'FontWeight', 'bold');
title('Halpha zoom', 'Color', 'k', 'FontWeight', 'bold');

nexttile([1 2]);
plot(ppm_axis, residual, '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 1.5);
hold on;
yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
style_axis(gca, 16);
xlim([1.1 4.1]);
ylim([-0.18 0.18]);
xlabel('^{1}H chemical shift (ppm)', 'Color', 'k', 'FontWeight', 'bold');
ylabel('fit - experiment', 'Color', 'k', 'FontWeight', 'bold');
title('Residual', 'Color', 'k', 'FontWeight', 'bold');

exportgraphics(fig, fullfile(project_dir, 'alanine_noesypr1d_experiment_fit.png'), ...
    'Resolution', 300, 'BackgroundColor', 'white');

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(project_dir, 'alanine_noesypr1d_experiment_fit.png'));
fprintf('  %s\n', fullfile(project_dir, 'alanine_noesypr1d_experiment_fit_summary.csv'));
fprintf('  %s\n', fullfile(project_dir, 'alanine_noesypr1d_experiment_fit_overlay.csv'));
fprintf('  %s\n', fullfile(project_dir, 'alanine_noesypr1d_experiment_fit.mat'));

%% Local helpers
function [rmse, details] = fit_rmse(q, fid_raw, SFO1_MHz, O1_Hz, SW_Hz, TD, ...
        ppm_axis, exp_ppm, exp_norm, score_mask, center_ppm)
    lb_Hz = exp(q(1));
    ppm_offset = q(2);
    phase_rad = q(3);

    if ~isfinite(lb_Hz) || lb_Hz < 0.01 || lb_Hz > 30 || ...
       ~isfinite(ppm_offset) || abs(ppm_offset) > 0.05 || ~isfinite(phase_rad)
        rmse = 1e6;
        details = struct();
        return
    end

    [~, spec_complex] = process_fid(fid_raw, SFO1_MHz, O1_Hz, SW_Hz, TD, lb_Hz);
    model_real = real(exp(1i*phase_rad) * spec_complex);
    target_full = interp1(exp_ppm + ppm_offset, exp_norm, ppm_axis, 'linear', 0);

    m = model_real(score_mask);
    target = target_full(score_mask);
    x = ppm_axis(score_mask) - center_ppm;

    good = isfinite(m) & isfinite(target);
    m = m(good);
    target = target(good);
    x = x(good);

    X = [m(:), ones(numel(m), 1), x(:)];
    coef = X \ target(:);
    fit = X * coef;

    rmse = sqrt(mean((fit - target(:)).^2));
    details = struct('coef', coef, 'target', target(:), 'fit', fit);
end

function [ppm_axis, spec_complex] = process_fid(fid, SFO1_MHz, O1_Hz, SW_Hz, TD, lb_Hz)
    fid = fid(:);
    dt = 1 / SW_Hz;
    t = (0:numel(fid)-1).' * dt;
    apod = exp(-pi * lb_Hz * t);
    fid_apod = fid .* apod;
    fid_apod(1) = fid_apod(1) / 2;
    spec_complex = fftshift(fft(fid_apod, TD));
    spec_complex = spec_complex(:);

    freq_Hz = linspace(-SW_Hz/2, SW_Hz/2, TD);
    ppm_axis = O1_Hz/SFO1_MHz - freq_Hz/SFO1_MHz;
    ppm_axis = fliplr(ppm_axis);
    ppm_axis = ppm_axis(:);
end

function [ppm_axis, y] = read_bruker_1r(path_1r)
    EXP_SF = 799.71;
    EXP_SW_p = 9615.38461538462;
    EXP_OFFSET = 10.71179;
    EXP_NC_proc = 10;
    EXP_BYTORDP = 0;
    EXP_DTYPP = 0;

    machinefmt = 'ieee-le';
    if EXP_BYTORDP ~= 0
        machinefmt = 'ieee-be';
    end

    fid = fopen(path_1r, 'r', machinefmt);
    if fid < 0
        error('Could not open %s', path_1r);
    end
    cleanup = onCleanup(@() fclose(fid));
    if EXP_DTYPP == 0
        y = fread(fid, inf, 'int32');
    else
        y = fread(fid, inf, 'double');
    end
    y = double(y(:)) * 2^EXP_NC_proc;
    clear cleanup;

    ppm_axis = EXP_OFFSET - (0:numel(y)-1).' * (EXP_SW_p / EXP_SF) / (numel(y)-1);
    [ppm_axis, order] = sort(ppm_axis);
    y = y(order);
end

function y = normalize_trace(y, norm_mask)
    y = y(:);
    y = y - median(y);
    sc = max(y(norm_mask));
    if ~isfinite(sc) || sc <= 0
        sc = max(abs(y));
    end
    if ~isfinite(sc) || sc <= 0
        sc = 1;
    end
    y = y / sc;
end

function mask = in_window(ppm_axis, win)
    mask = ppm_axis >= win(1) & ppm_axis <= win(2);
end

function r_value = trace_corr(a, b)
    a = a(:); b = b(:);
    a = a - mean(a);
    b = b - mean(b);
    denom = sqrt(sum(a.^2) * sum(b.^2));
    if denom <= eps
        r_value = NaN;
    else
        r_value = (a' * b) / denom;
    end
end

function deg = wrap_to_180(deg)
    deg = mod(deg + 180, 360) - 180;
end

function style_axis(ax, fs)
    set(ax, ...
        'Color', 'w', ...
        'XColor', 'k', ...
        'YColor', 'k', ...
        'GridColor', [0.83 0.83 0.83], ...
        'GridAlpha', 0.75, ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'LineWidth', 1.6, ...
        'Box', 'on', ...
        'XDir', 'reverse', ...
        'TickDir', 'out');
    grid(ax, 'on');
    xtickformat(ax, '%.2f');
end

function style_legend(lg)
    set(lg, ...
        'Color', 'w', ...
        'TextColor', 'k', ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.2, ...
        'FontSize', 15);
end
