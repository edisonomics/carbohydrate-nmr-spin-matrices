function results = simulate_alanine_spinach_fft(p)
% simulate_alanine_spinach_fft.m
%
% MULTI-FIELD GISSMO PIPELINE for the alanine AX3 spin system.
%
% Shared validated pipeline, factored out of the original per-field
% scripts (700 MHz BMRB, 600 MHz cece_data) so each new field is a
% short driver script, not a 500-line copy-paste. Do not duplicate this
% logic into a new script -- write a driver that calls this function.
%
% Validated decisions baked into this function (see field-specific
% driver scripts' history for how these were found):
%   - rho0/coil = state(spin_system,'L+','1H') (NOT hp_acquire's
%     Lz->Hermitian-Lx-pulse pattern, which produced a dispersive/
%     phase-distorted lineshape at both 700 and 600 MHz)
%   - ppm_axis = O1/SFO1 - freq_Hz/SFO1 (sign flip needed because the
%     rho0='L+' shortcut demodulates opposite to the standard NMR
%     receiver convention, mirroring the spectrum about the carrier)
%   - lb_Hz and a constant ppm calibration offset are jointly fit
%     against the EXPERIMENTAL trace (not the analytical theory curve)
%     -- a linewidth-only fit made r(Spinach,expt) WORSE at both fields
%     because the real problem was a small (<0.003 ppm) calibration
%     offset, not linewidth
%
% INPUT (struct p), required fields:
%   p.field_label     - e.g. '700 MHz' or '600 MHz' (used in titles/filenames)
%   p.project_dir     - working directory (cd'd into, outputs saved here)
%   p.spinach_root    - Spinach installation root
%   p.SFO1_MHz        - 1H observe frequency, MHz (acqus SFO1)
%   p.O1_Hz           - carrier offset from 0 ppm, Hz (acqus O1)
%   p.SW_Hz           - spectral width, Hz (acqus SW_h)
%   p.TD              - zero-fill target = experimental SI
%   p.exp_file        - full path to the processed Bruker '1r' file
%   p.EXP_SI          - procs SI
%   p.EXP_SF          - procs SF, MHz
%   p.EXP_SW_p        - procs SW_p, Hz
%   p.EXP_OFFSET      - procs OFFSET, ppm
%   p.EXP_BYTORDP     - procs BYTORDP (0 = little-endian)
%   p.EXP_DTYPP       - procs DTYPP (0 = int32, else double)
%   p.EXP_NC_proc     - procs NC_proc (scaling exponent)
%   p.output_suffix   - appended to all saved filenames, e.g. '700MHz'
%
% Optional fields (defaults shown):
%   p.delta_Ha_ppm    = 3.7680   (shared GISSMO matrix -- field-independent)
%   p.delta_CH3_ppm   = 1.4655
%   p.J_Ha_CH3_Hz     = 7.234
%   p.lb_Hz_guess     = 2.5      (only used as the fminsearch starting point)
%   p.methyl_win      = [1.30 1.65]
%   p.alpha_win       = [3.55 3.95]
%   p.exp_legend_label = p.field_label  (legend text for the experimental trace)
%
% OUTPUT: results struct with fitted lb_Hz, ppm_offset, all r/RMSE
% values, and the full model struct that gets saved to the .mat file.

%% ---- defaults ----
if ~isfield(p, 'delta_Ha_ppm'),    p.delta_Ha_ppm    = 3.7680;  end
if ~isfield(p, 'delta_CH3_ppm'),   p.delta_CH3_ppm   = 1.4655;  end
if ~isfield(p, 'J_Ha_CH3_Hz'),     p.J_Ha_CH3_Hz     = 7.234;   end
if ~isfield(p, 'lb_Hz_guess'),     p.lb_Hz_guess     = 2.5;     end
if ~isfield(p, 'methyl_win'),      p.methyl_win      = [1.30 1.65]; end
if ~isfield(p, 'alpha_win'),       p.alpha_win       = [3.55 3.95]; end
if ~isfield(p, 'exp_legend_label'), p.exp_legend_label = p.field_label; end

%% ============================================================
% 0. Paths
%% ============================================================

cd(p.project_dir);

if exist(p.spinach_root,'dir')
    % Put the requested Spinach checkout first if another copy is already
    % on the MATLAB path; this prevents an old create.m from being reused.
    addpath(genpath(p.spinach_root), '-begin');
end
workspace_root = fileparts(p.project_dir);
if exist(fullfile(workspace_root, 'spinach_pool_guard.m'), 'file') == 2
    addpath(workspace_root);
end
if exist('create','file') ~= 2
    error('Spinach create() not found. Check spinach_root.');
end
if ~isfield(p, 'parallel_workers') || isempty(p.parallel_workers)
    p.parallel_workers = 1;
else
    p.parallel_workers = max(1, round(p.parallel_workers));
end
spinach_pool_guard(p.parallel_workers, false);

%% ============================================================
% 1. Acquisition settings
%% ============================================================

SFO1_MHz = p.SFO1_MHz;
O1_Hz    = p.O1_Hz;
SW_Hz    = p.SW_Hz;
TD       = p.TD;
ncomplex = TD / 2;

B0_T = SFO1_MHz / 42.57747892;

fprintf('\n===== Acquisition settings (%s) =====\n', p.field_label);
fprintf('SFO1 = %.8f MHz  |  B0 = %.6f T\n', SFO1_MHz, B0_T);
fprintf('O1   = %.3f Hz  (= %.6f ppm)\n', O1_Hz, O1_Hz/SFO1_MHz);
fprintf('SW   = %.3f Hz  (= %.6f ppm)\n', SW_Hz,  SW_Hz/SFO1_MHz);

%% ============================================================
% 2. Spin system parameters (shared GISSMO matrix)
%% ============================================================

delta_Ha_ppm  = p.delta_Ha_ppm;
delta_CH3_ppm = p.delta_CH3_ppm;
J_Ha_CH3_Hz   = p.J_Ha_CH3_Hz;
lb_Hz_guess   = p.lb_Hz_guess;

fprintf('\n===== Spin system (shared GISSMO matrix) =====\n');
fprintf('Halpha shift = %.6f ppm\n', delta_Ha_ppm);
fprintf('CH3 shift    = %.6f ppm\n', delta_CH3_ppm);
fprintf('3J(HH)       = %.6f Hz\n', J_Ha_CH3_Hz);

%% ============================================================
% 3. Build Spinach spin system  (AX3: Halpha + 3 equivalent Hbeta)
%% ============================================================

sys.magnet   = B0_T;
sys.isotopes = {'1H','1H','1H','1H'};
sys.parallel = {'local', p.parallel_workers};

inter.zeeman.scalar = {delta_Ha_ppm, delta_CH3_ppm, delta_CH3_ppm, delta_CH3_ppm};

inter.coupling.scalar      = cell(4,4);
inter.coupling.scalar{1,2} = J_Ha_CH3_Hz;
inter.coupling.scalar{1,3} = J_Ha_CH3_Hz;
inter.coupling.scalar{1,4} = J_Ha_CH3_Hz;

bas.formalism     = 'sphten-liouv';
bas.approximation = 'none';

spin_system = create(sys, inter);
spin_system = basis(spin_system, bas);

%% ============================================================
% 4. Spinach pulse-acquire simulation -> FID
%% ============================================================

parameters.spins       = {'1H'};
parameters.offset      = O1_Hz;
parameters.sweep       = SW_Hz;
parameters.npoints     = ncomplex;
parameters.zerofill    = TD;
parameters.axis_units  = 'ppm';
parameters.invert_axis = 1;

% rho0='L+' shortcut -- validated baseline, see file header.
parameters.rho0 = state(spin_system, 'L+', '1H');
parameters.coil = state(spin_system, 'L+', '1H');

fprintf('\nRunning Spinach liquid() acquire...\n');
fid_raw = liquid(spin_system, @acquire, parameters, 'nmr');
fid_raw = fid_raw(:);

fprintf('FID length: %d complex points\n', numel(fid_raw));

dt = 1 / SW_Hz;

%% ============================================================
% 5. ppm axis construction (independent of lb_Hz)
%% ============================================================

freq_Hz  = linspace(-SW_Hz/2, SW_Hz/2, TD);
ppm_axis = O1_Hz/SFO1_MHz - freq_Hz/SFO1_MHz;     % sign flip, see header
ppm_axis = fliplr(ppm_axis);
ppm_axis = ppm_axis(:);

methyl_win = p.methyl_win;
alpha_win  = p.alpha_win;

ch3_mask_sp = ppm_axis >= methyl_win(1) & ppm_axis <= methyl_win(2);
ala_mask_native = (ppm_axis >= methyl_win(1) & ppm_axis <= methyl_win(2)) | ...
                  (ppm_axis >= alpha_win(1)  & ppm_axis <= alpha_win(2));

%% ============================================================
% 6. Read the processed 1r spectrum
%% ============================================================

if ~isfile(p.exp_file)
    error('Experimental 1r not found: %s', p.exp_file);
end

machinefmt = 'ieee-le';
if p.EXP_BYTORDP ~= 0, machinefmt = 'ieee-be'; end

fid_exp = fopen(p.exp_file, 'r', machinefmt);
if p.EXP_DTYPP == 0
    exp_y = fread(fid_exp, inf, 'int32');
else
    exp_y = fread(fid_exp, inf, 'double');
end
fclose(fid_exp);
exp_y = double(exp_y(:)) * 2^p.EXP_NC_proc;

fprintf('\nRead %d points from 1r  (expected %d)\n', numel(exp_y), p.EXP_SI);

exp_ppm = p.EXP_OFFSET - (0:numel(exp_y)-1).' * ...
          (p.EXP_SW_p / p.EXP_SF) / (numel(exp_y) - 1);

[exp_ppm_s, ix] = sort(exp_ppm);
exp_y_s = exp_y(ix);

exp_y_s = exp_y_s - median(exp_y_s);
ala_mask_exp = (exp_ppm_s >= methyl_win(1) & exp_ppm_s <= methyl_win(2)) | ...
               (exp_ppm_s >= alpha_win(1)  & exp_ppm_s <= alpha_win(2));

ch3_mask_exp = exp_ppm_s >= methyl_win(1) & exp_ppm_s <= methyl_win(2);
sc_exp = max(exp_y_s(ch3_mask_exp));
if sc_exp <= 0 || ~isfinite(sc_exp)
    sc_exp = max(abs(exp_y_s(ala_mask_exp)));
end
if sc_exp <= 0 || ~isfinite(sc_exp), sc_exp = 1; end

exp_norm = exp_y_s / sc_exp;

%% ============================================================
% 7. Jointly fit lb_Hz AND a constant ppm calibration offset
%    against the EXPERIMENTAL trace (not the theory curve)
%% ============================================================

joint_rmse_fn = @(q) joint_rmse_vs_expt( ...
    exp(q(1)), q(2), fid_raw, dt, TD, ch3_mask_sp, ala_mask_native, ...
    ppm_axis, exp_ppm_s, exp_norm);

q0 = [log(lb_Hz_guess), 0];
q_opt = fminsearch(joint_rmse_fn, q0);

lb_Hz      = exp(q_opt(1));
ppm_offset = q_opt(2);

fprintf('\n===== Joint linewidth + calibration fit (vs experiment) =====\n');
fprintf('Initial guess  : lb_Hz = %.3f Hz,  ppm_offset = 0\n', lb_Hz_guess);
fprintf('Fitted lb_Hz   : %.4f Hz\n', lb_Hz);
fprintf('Fitted offset  : %+.5f ppm  (added to the experimental axis)\n', ppm_offset);
fprintf('RMSE at fit    : %.5f\n', joint_rmse_fn(q_opt));

exp_ppm_s = exp_ppm_s + ppm_offset;

exp_on_sp = interp1(exp_ppm_s, exp_norm, ppm_axis, 'linear', 0);
ex_native = exp_on_sp(ala_mask_native);

%% ============================================================
% 8. Build the final Spinach spectrum at the fitted lb_Hz
%% ============================================================

[spec_spinach_norm, spec_spinach] = build_spinach_spectrum( ...
    fid_raw, dt, lb_Hz, TD, ch3_mask_sp);

fprintf('\n===== RAW SIGNAL DIAGNOSTIC =====\n');
fprintf('max(abs(fid_raw)) : %.6g\n', max(abs(fid_raw)));
[max_val, max_idx] = max(abs(spec_spinach));
fprintf('max(abs(spec_spinach)) : %.6g  at ppm %.4f  (index %d of %d)\n', ...
    max_val, ppm_axis(max_idx), max_idx, numel(ppm_axis));

tmp_ha  = ppm_axis(ppm_axis >= 3.60 & ppm_axis <= 3.90);
tmp_ch3 = ppm_axis(ppm_axis >= 1.35 & ppm_axis <= 1.60);
sp_ha   = spec_spinach(ppm_axis >= 3.60 & ppm_axis <= 3.90);
sp_ch3  = spec_spinach(ppm_axis >= 1.35 & ppm_axis <= 1.60);
if ~isempty(tmp_ha)
    [~,idx] = max(sp_ha);
    fprintf('Spinach Ha  max found at : %.4f ppm\n', tmp_ha(idx));
end
if ~isempty(tmp_ch3)
    [~,idx] = max(sp_ch3);
    fprintf('Spinach CH3 max found at : %.4f ppm\n', tmp_ch3(idx));
end

%% ============================================================
% 9. Build analytical AX3 first-order theory curve at the fitted lb_Hz
%% ============================================================

J_ppm    = J_Ha_CH3_Hz / SFO1_MHz;
fwhm_ppm = lb_Hz / SFO1_MHz;
hwhm_ppm = fwhm_ppm / 2;

lorentz = @(x, c, h) (h.^2) ./ ((x - c).^2 + h.^2);

ch3_lines = [delta_CH3_ppm - J_ppm/2; delta_CH3_ppm + J_ppm/2];
ha_lines  = [delta_Ha_ppm - 1.5*J_ppm; delta_Ha_ppm - 0.5*J_ppm; ...
             delta_Ha_ppm + 0.5*J_ppm; delta_Ha_ppm + 1.5*J_ppm];

ch3_weights = [1.5; 1.5];
ha_weights  = [1; 3; 3; 1] / 8;

spec_theory = zeros(size(ppm_axis));
for k = 1:2
    spec_theory = spec_theory + ch3_weights(k) * lorentz(ppm_axis, ch3_lines(k), hwhm_ppm);
end
for k = 1:4
    spec_theory = spec_theory + ha_weights(k) * lorentz(ppm_axis, ha_lines(k), hwhm_ppm);
end

sc_th = max(spec_theory(ch3_mask_sp));
if sc_th <= 0 || ~isfinite(sc_th), sc_th = 1; end

spec_theory_norm = spec_theory / sc_th;
spec_theory_norm = spec_theory_norm(:);

fprintf('\n===== Normalisation check =====\n');
fprintf('Spinach  CH3 peak after norm : %.4f  (should be 1.0)\n', max(spec_spinach_norm(ch3_mask_sp)));
fprintf('Theory   CH3 peak after norm : %.4f  (should be 1.0)\n', max(spec_theory_norm(ch3_mask_sp)));
fprintf('Expt     CH3 peak after norm : %.4f  (should be ~1.0)\n', max(exp_norm(ch3_mask_exp)));

%% ============================================================
% 10. Cross-validation diagnostics
%% ============================================================

sp_native = spec_spinach_norm(ala_mask_native);
th_native = spec_theory_norm(ala_mask_native);

good = isfinite(sp_native) & isfinite(th_native);
r_sp_th    = corr(sp_native(good), th_native(good));
rmse_sp_th = sqrt(mean((sp_native(good) - th_native(good)).^2));

good2     = isfinite(sp_native) & isfinite(ex_native);
r_sp_ex   = corr(sp_native(good2), ex_native(good2));
rmse_sp_ex = sqrt(mean((sp_native(good2) - ex_native(good2)).^2));

good3   = isfinite(th_native) & isfinite(ex_native);
r_th_ex = corr(th_native(good3), ex_native(good3));

fprintf('\n===== Cross-validation (native grid, fitted lb_Hz=%.4f Hz) =====\n', lb_Hz);
fprintf('Spinach vs theory   r = %.4f   RMSE = %.4f\n', r_sp_th, rmse_sp_th);
fprintf('Spinach vs expt     r = %.4f   RMSE = %.4f\n', r_sp_ex, rmse_sp_ex);
fprintf('Theory  vs expt     r = %.4f\n', r_th_ex);

%% ============================================================
% 11. Plots -- high-contrast, white background
%% ============================================================

close all;

col_exp     = [0.00 0.00 0.00];
col_spinach = [0.84 0.10 0.11];
col_theory  = [0.12 0.47 0.71];
col_marker  = [0.20 0.63 0.17];

suffix = p.output_suffix;

% ---- Fig 1: Full overlay ----
fig1 = figure('Color','w','Position',[60 60 1600 720]);
ax1 = axes('Parent', fig1);
hold(ax1, 'on');

h_exp = plot(ax1, exp_ppm_s, exp_norm, '-', ...
    'Color', col_exp, 'LineWidth', 2.5, ...
    'DisplayName', sprintf('%s experiment (%+.4f ppm calibrated)', p.exp_legend_label, ppm_offset));
h_sp = plot(ax1, ppm_axis, spec_spinach_norm, '--', ...
    'Color', col_spinach, 'LineWidth', 2.8, ...
    'DisplayName', sprintf('Spinach FFT  (LB = %.2f Hz, fitted)', lb_Hz));
h_th = plot(ax1, ppm_axis, spec_theory_norm, ':', ...
    'Color', col_theory, 'LineWidth', 3.0, ...
    'DisplayName', 'AX_3 analytical (first-order)');

style_ax(22);
xlim(ax1, [0.8 4.5]);
ylim(ax1, [-0.30 1.25]);
xlabel(ax1, '^{1}H chemical shift (ppm)', 'FontSize', 26, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax1, 'Normalised intensity',        'FontSize', 26, 'FontWeight', 'bold', 'Color', 'k');
title(ax1, { sprintf('Alanine %s  --  Experiment  |  Spinach FFT  |  AX_3 theory', p.field_label), ...
    sprintf('J = %.3f Hz    LB = %.3f Hz (fitted)    ppm offset = %+.4f (fitted)    r(theory) = %.4f    r(expt) = %.4f', ...
    J_Ha_CH3_Hz, lb_Hz, ppm_offset, r_sp_th, r_sp_ex) }, ...
    'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k');

legend(ax1, [h_exp h_sp h_th], 'Location', 'northwest', 'FontSize', 20, ...
    'TextColor', 'k', 'EdgeColor', 'k', 'Color', 'w', 'LineWidth', 1.2);

exportgraphics(fig1, sprintf('alanine_fft_overlay_full_%s.png', suffix), 'Resolution', 300);
fprintf('Saved alanine_fft_overlay_full_%s.png\n', suffix);

% ---- Fig 2: Zoom panels ----
fig2 = figure('Color','w','Position',[80 80 1600 700]);

ax2a = subplot(1, 2, 1);
hold(ax2a, 'on');
plot(ax2a, ppm_axis, spec_spinach_norm, '-', 'Color', col_spinach, 'LineWidth', 3.0, 'DisplayName', 'Spinach FFT');
plot(ax2a, ppm_axis, spec_theory_norm, '--', 'Color', col_theory, 'LineWidth', 3.0, 'DisplayName', 'AX_3 theory');
xline(ax2a, ch3_lines(1), '-', 'Color', col_marker, 'LineWidth', 2.0, 'Alpha', 0.9, ...
    'HandleVisibility', 'off', 'Label', sprintf('%.4f ppm', ch3_lines(1)), ...
    'LabelHorizontalAlignment', 'right', 'FontSize', 14, 'FontWeight', 'bold');
xline(ax2a, ch3_lines(2), '-', 'Color', col_marker, 'LineWidth', 2.0, 'Alpha', 0.9, ...
    'HandleVisibility', 'off', 'Label', sprintf('%.4f ppm', ch3_lines(2)), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 14, 'FontWeight', 'bold');
style_ax(20);
xlim(ax2a, methyl_win);
ylim(ax2a, [-0.20 1.20]);
xlabel(ax2a, '^{1}H chemical shift (ppm)', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax2a, 'Normalised intensity',        'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
title(ax2a,  'CH_3 doublet', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
legend(ax2a, 'Location', 'northwest', 'FontSize', 18, 'TextColor', 'k', 'EdgeColor', 'k', 'Color', 'w', 'LineWidth', 1.2);

ax2b = subplot(1, 2, 2);
hold(ax2b, 'on');
plot(ax2b, ppm_axis, spec_spinach_norm, '-', 'Color', col_spinach, 'LineWidth', 3.0, 'DisplayName', 'Spinach FFT');
plot(ax2b, ppm_axis, spec_theory_norm, '--', 'Color', col_theory, 'LineWidth', 3.0, 'DisplayName', 'AX_3 theory');
for k = 1:4
    xline(ax2b, ha_lines(k), '-', 'Color', col_marker, 'LineWidth', 2.0, 'Alpha', 0.9, ...
        'HandleVisibility', 'off', 'Label', sprintf('%.4f', ha_lines(k)), ...
        'LabelHorizontalAlignment', 'right', 'FontSize', 13, 'FontWeight', 'bold');
end
style_ax(20);
xlim(ax2b, alpha_win);
ylim(ax2b, [-0.10 0.55]);
xlabel(ax2b, '^{1}H chemical shift (ppm)', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax2b, 'Normalised intensity',        'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
title(ax2b,  'H_\alpha quartet', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
legend(ax2b, 'Location', 'northeast', 'FontSize', 18, 'TextColor', 'k', 'EdgeColor', 'k', 'Color', 'w', 'LineWidth', 1.2);

sgtitle(fig2, sprintf('%s: Spinach FFT vs AX_3 theory  |  r = %.4f  |  RMSE = %.4f  |  LB = %.3f Hz (fitted)', ...
    p.field_label, r_sp_th, rmse_sp_th, lb_Hz), 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');

exportgraphics(fig2, sprintf('alanine_fft_vs_theory_zoom_%s.png', suffix), 'Resolution', 300);
fprintf('Saved alanine_fft_vs_theory_zoom_%s.png\n', suffix);

% ---- Fig 3: Residual ----
fig3 = figure('Color','w','Position',[100 100 1600 480]);
ax3  = axes('Parent', fig3);
hold(ax3, 'on');

residual = spec_spinach_norm - spec_theory_norm;
ppm_row = ppm_axis(:)';
res_row = residual(:)';
fill(ax3, [ppm_row, fliplr(ppm_row)], [max(res_row, 0), zeros(1, numel(ppm_row))], ...
    col_spinach, 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
fill(ax3, [ppm_row, fliplr(ppm_row)], [min(res_row, 0), zeros(1, numel(ppm_row))], ...
    col_theory, 'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax3, ppm_axis, residual, '-', 'Color', 'k', 'LineWidth', 1.8, ...
    'DisplayName', 'Spinach FFT - AX_3 theory');
yline(ax3, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');

style_ax(20);
xlim(ax3, [0.8 4.5]);
ylim_max = max(0.05, max(abs(residual)) * 1.3);
ylim(ax3, [-ylim_max ylim_max]);
xlabel(ax3, '^{1}H chemical shift (ppm)', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax3, 'Residual (Spinach - theory)', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
title(ax3, sprintf('%s Residual  |  RMSE = %.5f  |  LB = %.3f Hz (fitted)', ...
    p.field_label, rmse_sp_th, lb_Hz), 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
legend(ax3, 'Location', 'northwest', 'FontSize', 18, 'TextColor', 'k', 'EdgeColor', 'k', 'Color', 'w', 'LineWidth', 1.2);

exportgraphics(fig3, sprintf('alanine_fft_residual_%s.png', suffix), 'Resolution', 300);
fprintf('Saved alanine_fft_residual_%s.png\n', suffix);

%% ============================================================
% 12. Save outputs
%% ============================================================

model.field_label  = p.field_label;
model.spin_names   = {'Halpha','Hbeta1','Hbeta2','Hbeta3'};
model.shift_ppm    = [delta_Ha_ppm; delta_CH3_ppm; delta_CH3_ppm; delta_CH3_ppm];
model.J_Hz         = J_Ha_CH3_Hz;
model.lb_Hz        = lb_Hz;
model.lb_Hz_guess  = lb_Hz_guess;
model.ppm_offset_fitted = ppm_offset;
model.SFO1_MHz     = SFO1_MHz;
model.O1_Hz        = O1_Hz;
model.SW_Hz        = SW_Hz;
model.B0_T         = B0_T;
model.r_spinach_vs_theory = r_sp_th;
model.r_spinach_vs_expt   = r_sp_ex;
model.rmse_spinach_vs_theory = rmse_sp_th;
model.rmse_spinach_vs_expt   = rmse_sp_ex;
model.r_theory_vs_expt    = r_th_ex;

save(sprintf('alanine_spinach_fft_model_%s.mat', suffix), 'model');

overlay = table(ppm_axis(:), spec_spinach_norm(:), spec_theory_norm(:), exp_on_sp(:), ...
    'VariableNames', {'ppm','spinach_norm','theory_norm','experiment_interp_norm'});
writetable(overlay, sprintf('alanine_spinach_fft_overlay_%s.csv', suffix));

fprintf('\nSaved:\n');
fprintf('  alanine_fft_overlay_full_%s.png\n', suffix);
fprintf('  alanine_fft_vs_theory_zoom_%s.png\n', suffix);
fprintf('  alanine_fft_residual_%s.png\n', suffix);
fprintf('  alanine_spinach_fft_model_%s.mat\n', suffix);
fprintf('  alanine_spinach_fft_overlay_%s.csv\n', suffix);

results = model;

end

%% ============================================================
% Local helper functions
%% ============================================================

% ---- Apodize + FFT + CH3-normalize a Spinach FID at a given lb_Hz ----
function [spec_norm, spec_unnorm] = build_spinach_spectrum( ...
    fid_raw, dt, lb_Hz, TD, ch3_mask_sp)

    t    = (0:numel(fid_raw)-1).' * dt;
    apod = exp(-pi * lb_Hz * t);
    fid_apod = fid_raw .* apod;
    fid_apod(1) = fid_apod(1) / 2;

    spec_raw = fftshift(fft(fid_apod, TD));
    spec = real(spec_raw);
    spec = spec(:);
    % NOTE: do NOT fliplr() here -- spec is a column vector, and fliplr
    % on a column is a documented no-op. ppm_axis IS flipped (as a row,
    % before being columnized) in section 5, and that asymmetry nets
    % out to the correct pairing empirically (r > 0.99 at two fields).
    % See the original 700/600 MHz scripts' version history.

    bl = median(spec);
    spec_bl = spec - bl;

    sc = max(spec_bl(ch3_mask_sp));
    if sc <= 0 || ~isfinite(sc), sc = 1; end

    spec_norm   = spec_bl / sc;
    spec_unnorm = spec;
end

% ---- RMSE of Spinach spectrum vs experiment, jointly in lb_Hz and a
% ---- constant ppm offset applied to the experimental axis.
function r = joint_rmse_vs_expt( ...
    lb_Hz, ppm_offset, fid_raw, dt, TD, ch3_mask_sp, ala_mask_native, ...
    ppm_axis, exp_ppm_s, exp_norm)

    exp_on_sp_trial = interp1(exp_ppm_s + ppm_offset, exp_norm, ppm_axis, ...
        'linear', 0);
    ex_native_trial = exp_on_sp_trial(ala_mask_native);

    spec_norm = build_spinach_spectrum(fid_raw, dt, lb_Hz, TD, ch3_mask_sp);
    sp_native = spec_norm(ala_mask_native);

    r = sqrt(mean((sp_native - ex_native_trial).^2));
end

% ---- Shared axis-formatting helper -------------------------
function style_ax(fs)
    ax = gca;
    set(ax, ...
        'Color',           'w', ...
        'XColor',          'k', ...
        'YColor',          'k', ...
        'GridColor',       [0.88 0.88 0.88], ...
        'GridAlpha',       0.6, ...
        'MinorGridColor',  [0.88 0.88 0.88], ...
        'FontSize',        fs, ...
        'FontWeight',      'bold', ...
        'LineWidth',       1.8, ...
        'Box',             'on', ...
        'XDir',            'reverse', ...
        'TickDir',         'out', ...
        'TickLength',      [0.012 0.012]);
    grid on;
    ax.Toolbar.Visible = 'off';
end
