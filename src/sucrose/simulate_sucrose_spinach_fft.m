function results = simulate_sucrose_spinach_fft(p)
% simulate_sucrose_spinach_fft.m
%
% N-SPIN GENERALIZATION of the alanine pipeline (simulate_alanine_spinach_fft.m)
% for sucrose's 14-proton GISSMO matrix (bmse000119, glucose+fructose rings,
% strongly coupled -- no first-order analytical theory equivalent exists for
% this system, see simulate_alanine_spinach_fft.m header for why that
% matters).
%
% Carries forward from the alanine work:
%   - direct mode uses rho0/coil = state(spin_system,'L+','1H') as a
%     matrix-only shortcut (validated, NOT hp_acquire)
%   - noesypr1d mode uses equilibrium rho0=Lz and L+ detection coil
%   - ppm_axis = O1/SFO1 - freq_Hz/SFO1, fliplr ppm only, NOT spec (see
%     build_sucrose_spectrum below for the no-op-fliplr-on-column note)
%
% NEW relative to the alanine function (sucrose-specific generalizations):
%   - N spins read directly from a GISSMO matrix file (not hardcoded AX3)
%   - bas.approximation = 'IK-2' (full Liouville space intractable at 14
%     spins), bas.connectivity = 'scalar_couplings', bas.space_level = 3
%     -- these three settings are carried over from a prior independent
%     attempt at this same system (simulate_sucrose_from_GISSMO_matrix.m,
%     see sucrose_refinement_notes.txt), which got the simulation running
%     at this spin count with them.
%   - ncomplex (acquired complex points) and TD (zero-fill target) are
%     INDEPENDENT parameters, not assumed TD/2 -- this dataset acquired
%     32768 complex points and zero-filled to 32768 (ratio 1), unlike
%     both alanine datasets (ratio 2).
%
% DELIBERATELY NOT carried over from the prior sucrose attempt:
%   - Its ppm-axis fix (offset_ppm + hz_axis/MHz, no fliplr) -- ad hoc for
%     its own parameterization, not connected to the fftshift/freq_Hz
%     construction validated twice on alanine (r>0.97 at two fields).
%   - Its hand-edited proton shifts (e.g. moving one shift from 3.734 to
%     5.221 ppm to chase a peak match) -- conflates "is the pipeline
%     right" with "is the matrix right". This function uses the
%     UNEDITED official GISSMO matrix only.
%
% INPUT (struct p), required fields:
%   p.matrix_file   - path to the 14x14 numeric GISSMO matrix (diagonal =
%                     shifts in ppm, off-diagonal = J couplings in Hz)
%   p.atom_ids      - cell array of the 14 GISSMO atom-ID labels, in the
%                     same row/column order as the matrix file (for
%                     traceability back to the bmse000119 entry only --
%                     not used in any physics)
%   p.project_dir, p.spinach_root
%   p.SFO1_MHz, p.O1_Hz, p.SW_Hz
%   p.ncomplex      - acquired complex points (acqus TD / 2)
%   p.TD            - zero-fill target (procs SI)
%   p.lb_Hz_guess
%
% Optional:
%   p.blocks            - cell array of 1-based index vectors into the
%                          14-spin matrix, each defining an INDEPENDENT
%                          (zero cross-coupling) block, e.g.
%                          {[1:7],[8:14]} for sucrose's glucose/fructose
%                          rings. When given, each block is simulated
%                          EXACTLY (bas.approximation='none', tractable
%                          since each block is small) and the FIDs are
%                          summed (exact, not approximate, because the
%                          cross-coupling truly is zero). When omitted,
%                          falls back to a single IK-2-approximated run
%                          over all spins together.
%   p.bas_approximation (default 'IK-2', ignored if p.blocks is given)
%   p.bas_space_level   (default 3, ignored if p.blocks is given)
%   p.anomeric_win       (default [5.30 5.50]) -- normalization/sanity window
%   p.water_win          (default [4.65 4.90]) -- excluded from the
%                          experiment FIT only (HOD/water peak)
%   p.artifact_win       (default [5.15 5.30]) -- excluded from the
%                          experiment FIT only (reproducible unexplained
%                          peak documented in sucrose_refinement_notes.txt,
%                          not predicted by GISSMO or prominent in BMRB)
%   p.crowded_win        (default [3.60 3.73]) -- NOT excluded from the
%                          fit, only from the with/without diagnostic
%                          comparison reported after fitting. Brackets
%                          the matrix's near-degenerate 3.667/3.670 ppm
%                          shifts, where a small per-sample shift
%                          difference shows up as a large localized
%                          residual even when the rest of the spectrum
%                          aligns well.
%   p.gissmo_sim_file    - path to GISSMO's published sim_<field>MHz.json
%                          (e.g. .../spectral_data/sim_1100MHz.json). The
%                          "theory" tier for sucrose -- see header notes.
%   p.exp_file           - path to a processed Bruker '1r' file (real
%                          experimental spectrum). The AUTHORITATIVE fit
%                          target when given: lb_Hz and a ppm calibration
%                          offset are fit against this, not GISSMO's
%                          simulation, mirroring how alanine fit against
%                          real BMRB/lab data with theory recomputed
%                          afterward at the same linewidth. Needs
%                          p.EXP_SI, EXP_SF, EXP_SW_p, EXP_OFFSET,
%                          EXP_BYTORDP, EXP_DTYPP, EXP_NC_proc (same
%                          meaning as the alanine scripts' procs fields).
%                          If omitted, falls back to fitting against
%                          GISSMO's simulation instead.
%   p.parallel_workers   - MATLAB workers for Spinach create() (default 1).
%                          Stale pools are cleaned once per MATLAB session
%                          by spinach_pool_guard before the first create().
%   p.suppress_sucrose_diagnostics - suppress hard-coded sucrose labels in
%                          the raw diagnostic (default false).
%   p.diagnostic_peaks   - optional struct array with fields label and win;
%                          replaces the legacy sucrose peak checks.
%
% OUTPUT: results struct with the fitted lb_Hz, both ppm calibration
% offsets (vs experiment, vs GISSMO), and the three-way r/RMSE
% (Spinach-vs-experiment, Spinach-vs-GISSMO, GISSMO-vs-experiment).

% Defaults come from the checked-in molecule configuration. Callers may
% override fields in p, but the workflow settings are not hidden literals.
src_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(src_dir));
if ~isfield(p, 'molecule') || isempty(p.molecule), p.molecule = 'sucrose'; end
if exist(fullfile(repo_dir, 'src', 'common', 'load_carbohydrate_config.m'), 'file')
    addpath(fullfile(repo_dir, 'src', 'common'));
    cfg = load_carbohydrate_config(repo_dir, p.molecule);
else
    cfg = struct();
end
if isfield(cfg, 'spinach')
    if ~isfield(p, 'bas_approximation'), p.bas_approximation = cfg.spinach.basis_approximation; end
    if ~isfield(p, 'bas_connectivity'), p.bas_connectivity = cfg.spinach.basis_connectivity; end
    if ~isfield(p, 'bas_space_level'), p.bas_space_level = cfg.spinach.basis_space_level; end
end
if isfield(cfg, 'sequence')
    if ~isfield(p, 'parallel_workers'), p.parallel_workers = cfg.sequence.parallel_workers; end
    if ~isfield(p, 'noesy_tmix_s'), p.noesy_tmix_s = cfg.sequence.tmix_s; end
    if ~isfield(p, 'noesy_d1_s'), p.noesy_d1_s = cfg.sequence.d1_s; end
end
if isfield(cfg, 'processing')
    if ~isfield(p, 'water_win'), p.water_win = cfg.processing.water_region_ppm; end
    if ~isfield(p, 'artifact_win'), p.artifact_win = cfg.processing.artifact_region_ppm; end
    if ~isfield(p, 'anomeric_win'), p.anomeric_win = cfg.processing.anomeric_region_ppm; end
    if ~isfield(p, 'crowded_win'), p.crowded_win = cfg.processing.crowded_region_ppm; end
    if ~isfield(p, 'sucrose_region'), p.sucrose_region = cfg.processing.fit_region_ppm; end
end

if ~isfield(p, 'bas_approximation'), p.bas_approximation = 'IK-2'; end
if ~isfield(p, 'bas_space_level'),   p.bas_space_level   = 3;      end
if ~isfield(p, 'sample_label') || isempty(p.sample_label)
    p.sample_label = p.molecule;
end
if ~isfield(p, 'field_label') || isempty(p.field_label)
    p.field_label = 'field';
end
if ~isfield(p, 'plot_prefix') || isempty(p.plot_prefix)
    p.plot_prefix = [p.molecule '_spinach'];
end
% Spinach's create() starts/reuses a parallel pool by default.  Use one
% worker unless a caller explicitly requests another count; the shared guard
% clears stale state once per MATLAB session before the first create().
if ~isfield(p, 'parallel_workers') || isempty(p.parallel_workers)
    p.parallel_workers = 1;
else
    p.parallel_workers = max(1, round(p.parallel_workers));
end
if ~isfield(p, 'acquisition_mode') || isempty(p.acquisition_mode)
    p.acquisition_mode = 'direct';
end
if ~isfield(p, 'noesy_tmix_s') || isempty(p.noesy_tmix_s)
    p.noesy_tmix_s = 0.05;
end
if ~isfield(p, 'noesy_d1_s') || isempty(p.noesy_d1_s)
    p.noesy_d1_s = 2.0;
end
if ~isfield(p, 'noesy_presat_nu1_Hz') || isempty(p.noesy_presat_nu1_Hz)
    p.noesy_presat_nu1_Hz = 0;
end
if ~isfield(p, 'noesy_presat_d1') || isempty(p.noesy_presat_d1)
    p.noesy_presat_d1 = false;
end
if ~isfield(p, 'noesy_presat_tmix') || isempty(p.noesy_presat_tmix)
    p.noesy_presat_tmix = false;
end
if ~isfield(p, 'noesy_pulse_sign') || isempty(p.noesy_pulse_sign)
    p.noesy_pulse_sign = 1;
end
if ~isfield(p, 'noesy_receiver_sign') || isempty(p.noesy_receiver_sign)
    p.noesy_receiver_sign = 1;
end
if ~isfield(p, 'noesy_phase_cycle') || isempty(p.noesy_phase_cycle)
    p.noesy_phase_cycle = true;
end
if ~isfield(p, 'noesy_sequence_model') || isempty(p.noesy_sequence_model)
    p.noesy_sequence_model = 'homospoil';
end
if ~isfield(p, 'noesy_storage_filter') || isempty(p.noesy_storage_filter)
    p.noesy_storage_filter = 'zeroq';
end
if ~isfield(p, 'noesy_crusher_mode') || isempty(p.noesy_crusher_mode)
    p.noesy_crusher_mode = 'none';
end
if ~isfield(p, 'noesy_final_pulse_angle') || isempty(p.noesy_final_pulse_angle)
    p.noesy_final_pulse_angle = pi/2;
end
acq_mode = lower(strtrim(p.acquisition_mode));
if ~ismember(acq_mode, {'direct', 'noesypr1d'})
    error('Unknown acquisition_mode "%s"; use direct or noesypr1d.', p.acquisition_mode);
end
if ~isfield(p, 'fit_receiver_phase') || isempty(p.fit_receiver_phase)
    p.fit_receiver_phase = strcmp(acq_mode, 'noesypr1d');
end
fit_receiver_phase = logical(p.fit_receiver_phase);
if ~isfield(p, 'lineshape_model') || isempty(p.lineshape_model)
    p.lineshape_model = 'voigt';
end
lineshape_model = lower(strtrim(p.lineshape_model));
if ~ismember(lineshape_model, {'voigt', 'lorentzian'})
    error('Unknown lineshape_model "%s"; use voigt or lorentzian.', p.lineshape_model);
end

%% ---- Paths ----
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
spinach_pool_guard(p.parallel_workers, false);

%% ---- Load the GISSMO matrix ----
M = readmatrix(p.matrix_file, 'FileType', 'text');
if size(M,1) ~= size(M,2)
    error('GISSMO matrix is not square (got %dx%d).', size(M,1), size(M,2));
end
nspins = size(M,1);

shifts_ppm = diag(M);
J_Hz = triu(M,1);
J_Hz = J_Hz + J_Hz.';   % symmetrize (file stores upper triangle only)

% ---- Optional per-spin INTRINSIC linewidth (two-tier model) ----
% p.lb_intrinsic_Hz, if given, is an nspins-vector of intrinsic FWHM (Hz)
% per spin (matrix/atom_ids order), applied INSIDE the Spinach simulation
% as a phenomenological t1_t2 relaxation superoperator (R2 = pi*FWHM per
% spin -- the exact equivalence to apodization proved in
% demo_linewidth_mechanisms_alanine_600MHz.m). Physically these are the
% dipolar relaxation widths from demo_redfield_intrinsic_linewidth_sucrose.m
% (CH2/CH2OH protons ~0.75 Hz, ring CH ~0.18 Hz). When present, the flat
% lb_Hz fit downstream then only carries the FIELD/SHIM INHOMOGENEITY that
% is common to all protons, so total width per proton = intrinsic(k) + lb.
% When ABSENT, behaviour is identical to the validated flat-lb baseline.
have_intrinsic = isfield(p,'lb_intrinsic_Hz') && ~isempty(p.lb_intrinsic_Hz);
if have_intrinsic
    lb_intrinsic = p.lb_intrinsic_Hz(:).';
    if numel(lb_intrinsic) ~= nspins
        error('p.lb_intrinsic_Hz has %d entries but the matrix has %d spins.', ...
            numel(lb_intrinsic), nspins);
    end
    fprintf('\nTWO-TIER intrinsic linewidths ON: %d values, range %.3f-%.3f Hz.\n', ...
        nspins, min(lb_intrinsic), max(lb_intrinsic));
end

fprintf('\n===== GISSMO matrix =====\n');
fprintf('Loaded %d-spin matrix from %s\n', nspins, p.matrix_file);
if isfield(p, 'atom_ids') && numel(p.atom_ids) == nspins
    for k = 1:nspins
        fprintf('  spin %2d (atom %s): shift = %.4f ppm\n', k, p.atom_ids{k}, shifts_ppm(k));
    end
end

%% ---- Acquisition settings ----
SFO1_MHz = p.SFO1_MHz;
O1_Hz    = p.O1_Hz;
SW_Hz    = p.SW_Hz;
ncomplex = p.ncomplex;
TD       = p.TD;

gamma_mhz_per_t = 42.57747892;
if isfield(cfg, 'spinach') && isfield(cfg.spinach, 'proton_gamma_mhz_per_t')
    gamma_mhz_per_t = cfg.spinach.proton_gamma_mhz_per_t;
end
B0_T = SFO1_MHz / gamma_mhz_per_t;

fprintf('\n===== Acquisition settings =====\n');
fprintf('SFO1 = %.8f MHz  |  B0 = %.6f T\n', SFO1_MHz, B0_T);
fprintf('O1   = %.3f Hz  (= %.6f ppm)\n', O1_Hz, O1_Hz/SFO1_MHz);
fprintf('SW   = %.3f Hz  (= %.6f ppm)\n', SW_Hz, SW_Hz/SFO1_MHz);
fprintf('ncomplex (acquired) = %d,  TD (zerofill target) = %d\n', ncomplex, TD);
fprintf('Acquisition mode = %s', acq_mode);
if strcmp(acq_mode, 'noesypr1d')
    fprintf(['  (d1 %.3f s, tmix %.3f s, pulse sign %+g, receiver sign %+g, ' ...
        'model %s, phase cycle %d, storage %s, crusher %s, presat %.2f Hz d1=%d tmix=%d)\n'], ...
        p.noesy_d1_s, p.noesy_tmix_s, p.noesy_pulse_sign, p.noesy_receiver_sign, ...
        p.noesy_sequence_model, logical(p.noesy_phase_cycle), ...
        p.noesy_storage_filter, p.noesy_crusher_mode, ...
        p.noesy_presat_nu1_Hz, ...
        logical(p.noesy_presat_d1), logical(p.noesy_presat_tmix));
    fprintf('  Initial state = Lz; detection coil = L+ (L+ is not used as rho0).\n');
else
    fprintf('\n');
end
fprintf('Line-shape fit model = %s\n', lineshape_model);

%% ---- Build and run the Spinach spin system(s) ----
% Two modes:
%   p.blocks given: the matrix decomposes into independent blocks with
%     ZERO cross-coupling (true for sucrose: glucose-ring spins 1-7 and
%     fructose-ring spins 8-14 never couple, since the glycosidic link
%     is a 5+-bond H...H pathway -- see conversation notes). Each block
%     is small enough to simulate EXACTLY (bas.approximation='none', no
%     truncation at all), and because cross-coupling is exactly zero,
%     the joint FID is exactly the SUM of each block's independently-
%     simulated FID (the evolution operator factorizes as a tensor
%     product across blocks, and the L+ coil sums linearly over them --
%     no approximation in this decomposition itself, only in dropping
%     the (zero) cross-terms, which is exact, not approximate).
%   p.blocks omitted: single IK-2-approximated run over all spins
%     together (the original approach, kept for molecules that don't
%     decompose this cleanly).

if isfield(p, 'blocks') && ~isempty(p.blocks)

    fprintf('\nRunning Spinach as %d independent EXACT block(s) (bas=''none''):\n', ...
        numel(p.blocks));
    fid_raw = zeros(ncomplex, 1);
    tic;
    for bk = 1:numel(p.blocks)
        idx = p.blocks{bk};
        nb  = numel(idx);

        sys_b.magnet   = B0_T;
        sys_b.isotopes = repmat({'1H'}, 1, nb);
        if isfield(p, 'parallel_workers') && ~isempty(p.parallel_workers)
            sys_b.parallel = {'local', p.parallel_workers};
        end

        inter_b.zeeman.scalar = cell(1, nb);
        for k = 1:nb
            inter_b.zeeman.scalar{k} = shifts_ppm(idx(k));
        end

        % Spinach scalar J entries are supplied once per pair.
        inter_b.coupling.scalar = cell(nb, nb);
        for a = 1:nb-1
            for b = a+1:nb
                if J_Hz(idx(a), idx(b)) ~= 0
                    inter_b.coupling.scalar{a,b} = J_Hz(idx(a), idx(b));
                end
            end
        end

        % Per-spin intrinsic linewidth (two-tier) as an in-sim relaxation
        % superoperator -- same t1_t2 settings validated in the alanine
        % linewidth demo. R1 is set equal to R2 (harmless: only the
        % transverse R2 sets the detected L+ lineshape).
        if have_intrinsic
            r2b = pi * lb_intrinsic(idx);            % s^-1, per spin in block
            inter_b.relaxation  = {'t1_t2'};
            inter_b.equilibrium = 'zero';
            inter_b.rlx_keep    = 'diagonal';
            inter_b.temperature = 298;
            inter_b.r1_rates    = num2cell(r2b);
            inter_b.r2_rates    = num2cell(r2b);
        end

        bas_b.formalism     = 'sphten-liouv';
        bas_b.approximation = 'none';   % exact -- tractable at this block size

        spin_system_b = create(sys_b, inter_b);
        spin_system_b = basis(spin_system_b, bas_b);

        parameters_b.spins       = {'1H'};
        parameters_b.offset      = O1_Hz;
        parameters_b.sweep       = SW_Hz;
        parameters_b.npoints     = ncomplex;
        parameters_b.zerofill    = TD;
        parameters_b.axis_units  = 'ppm';
        parameters_b.invert_axis = 1;
        parameters_b.decouple    = {};
        parameters_b.coil = state(spin_system_b, 'L+', '1H');
        if strcmp(acq_mode, 'noesypr1d')
            parameters_b.rho0 = state(spin_system_b, 'Lz', '1H', 'cheap');
            parameters_b.tmix = p.noesy_tmix_s;
            parameters_b.d1 = p.noesy_d1_s;
            parameters_b.presat_nu1_Hz = p.noesy_presat_nu1_Hz;
            parameters_b.presat_d1 = p.noesy_presat_d1;
            parameters_b.presat_tmix = p.noesy_presat_tmix;
            parameters_b.pulse_sign = p.noesy_pulse_sign;
            parameters_b.receiver_sign = p.noesy_receiver_sign;
            parameters_b.phase_cycle = p.noesy_phase_cycle;
            parameters_b.sequence_model = p.noesy_sequence_model;
            parameters_b.storage_filter = p.noesy_storage_filter;
            parameters_b.crusher_mode = p.noesy_crusher_mode;
            parameters_b.final_pulse_angle = p.noesy_final_pulse_angle;
            experiment_fn = @noesypr1d_acquire;
        else
            parameters_b.rho0 = state(spin_system_b, 'L+', '1H');
            experiment_fn = @acquire;
        end

        fprintf('  block %d: %d spins (atoms', bk, nb);
        if isfield(p, 'atom_ids')
            fprintf(' %s', strjoin(p.atom_ids(idx), ','));
        end
        fprintf(')...\n');

        fid_b = liquid(spin_system_b, experiment_fn, parameters_b, 'nmr');
        fid_raw = fid_raw + fid_b(:);
    end
    fprintf('Done in %.1f s (all blocks). FID length: %d complex points\n', toc, numel(fid_raw));

else
    sys.magnet   = B0_T;
    sys.isotopes = repmat({'1H'}, 1, nspins);
    if isfield(p, 'parallel_workers') && ~isempty(p.parallel_workers)
        sys.parallel = {'local', p.parallel_workers};
    end

    inter.zeeman.scalar = cell(1, nspins);
    for k = 1:nspins
        inter.zeeman.scalar{k} = shifts_ppm(k);
    end

    % Spinach scalar J entries are supplied once per pair.
    inter.coupling.scalar = cell(nspins, nspins);
    for a = 1:nspins-1
        for b = a+1:nspins
            if J_Hz(a,b) ~= 0
                inter.coupling.scalar{a,b} = J_Hz(a,b);
            end
        end
    end

    if have_intrinsic
        r2a = pi * lb_intrinsic;                 % s^-1, per spin
        inter.relaxation  = {'t1_t2'};
        inter.equilibrium = 'zero';
        inter.rlx_keep    = 'diagonal';
        inter.temperature = 298;
        inter.r1_rates    = num2cell(r2a);
        inter.r2_rates    = num2cell(r2a);
    end

    bas.formalism     = 'sphten-liouv';
    bas.approximation = p.bas_approximation;
    if isfield(p, 'bas_connectivity')
        bas.connectivity = p.bas_connectivity;
    else
        bas.connectivity = 'scalar_couplings';
    end
    bas.space_level    = p.bas_space_level;

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
    parameters.coil = state(spin_system, 'L+', '1H');

    if strcmp(acq_mode, 'noesypr1d')
        parameters.rho0 = state(spin_system, 'Lz', '1H', 'cheap');
        parameters.tmix = p.noesy_tmix_s;
        parameters.d1 = p.noesy_d1_s;
        parameters.presat_nu1_Hz = p.noesy_presat_nu1_Hz;
        parameters.presat_d1 = p.noesy_presat_d1;
        parameters.presat_tmix = p.noesy_presat_tmix;
        parameters.pulse_sign = p.noesy_pulse_sign;
        parameters.receiver_sign = p.noesy_receiver_sign;
        parameters.phase_cycle = p.noesy_phase_cycle;
        parameters.sequence_model = p.noesy_sequence_model;
        parameters.storage_filter = p.noesy_storage_filter;
        parameters.crusher_mode = p.noesy_crusher_mode;
        parameters.final_pulse_angle = p.noesy_final_pulse_angle;
        experiment_fn = @noesypr1d_acquire;
    else
        % rho0='L+' shortcut -- validated on alanine at two fields, NOT hp_acquire.
        parameters.rho0 = state(spin_system, 'L+', '1H');
        experiment_fn = @acquire;
    end

    fprintf('\nRunning Spinach liquid() %s (%d spins, basis=%s)...\n', ...
        acq_mode, nspins, p.bas_approximation);
    tic;
    fid_raw = liquid(spin_system, experiment_fn, parameters, 'nmr');
    fid_raw = fid_raw(:);
    fprintf('Done in %.1f s. FID length: %d complex points\n', toc, numel(fid_raw));
end

% Label for plots/legends -- reflects which method actually ran.
if isfield(p, 'blocks') && ~isempty(p.blocks)
    method_label = sprintf('exact, %d blocks', numel(p.blocks));
else
    method_label = p.bas_approximation;
end

dt = 1 / SW_Hz;

%% ---- ppm axis construction (validated alanine formula, unchanged) ----
freq_Hz  = linspace(-SW_Hz/2, SW_Hz/2, TD);
ppm_axis = O1_Hz/SFO1_MHz - freq_Hz/SFO1_MHz;
ppm_axis = fliplr(ppm_axis);
ppm_axis = ppm_axis(:);

%% ---- Anomeric anchor window (normalization + sanity check) ----
% Glucose anomeric H1 (atom 37, shift 5.403 ppm) is large, isolated, and
% away from both the 4.65-4.90 water band and the 5.15-5.30 unexplained-
% artifact band documented in sucrose_refinement_notes.txt. Used here
% both as the axis-sign sanity check AND as the normalization anchor
% (replacing global-max normalization, per the prior attempt's lesson
% that global-max is fragile if an unexpected artifact happens to be
% tallest).
if isfield(p, 'anomeric_win'), anomeric_win = p.anomeric_win; else, anomeric_win = [5.30 5.50]; end
anom_mask = ppm_axis >= anomeric_win(1) & ppm_axis <= anomeric_win(2);
if ~any(anom_mask)
    error('No ppm points fell inside the anomeric anchor window [%.2f %.2f].', ...
        anomeric_win(1), anomeric_win(2));
end

fprintf('\n===== RAW SIGNAL DIAGNOSTIC =====\n');
fprintf('max(abs(fid_raw))      : %.6g\n', max(abs(fid_raw)));
[~, spec_check] = build_sucrose_spectrum(fid_raw, dt, p.lb_Hz_guess, 0, 0, TD, anom_mask);
[max_val, max_idx] = max(abs(spec_check));
fprintf('max(abs(spec_spinach)) : %.6g  at ppm %.4f  (index %d of %d)\n', ...
    max_val, ppm_axis(max_idx), max_idx, numel(ppm_axis));
tmp_ppm = ppm_axis(anom_mask);
tmp_y   = spec_check(anom_mask);
[~, idx] = max(abs(tmp_y));
if isfield(p, 'suppress_sucrose_diagnostics') && p.suppress_sucrose_diagnostics
    fprintf('Custom-matrix anchor max found at : %.4f ppm\n', tmp_ppm(idx));
else
    fprintf('Glucose anomeric H1 max found at : %.4f ppm  (expected 5.403 ppm)\n', tmp_ppm(idx));
end

%% ============================================================
% Comparison regions and exclusion windows
%% ============================================================
% sucrose_region: overall comparison window, per sucrose_refinement_notes.txt.
% water_win/artifact_win: excluded from fitting against REAL experimental
% data only (GISSMO's simulation has neither feature, so its comparison
% uses the full sucrose_region with no exclusions).

sucrose_region = [3.0 5.8];
if isfield(p, 'sucrose_region'), sucrose_region = p.sucrose_region; end
if isfield(p, 'water_win'),    water_win    = p.water_win;    else, water_win    = [4.65 4.90]; end
if isfield(p, 'artifact_win'), artifact_win = p.artifact_win; else, artifact_win = [5.15 5.30]; end

region_mask = ppm_axis >= sucrose_region(1) & ppm_axis <= sucrose_region(2);
fit_mask_expt = region_mask;
if numel(water_win) >= 2
    fit_mask_expt = fit_mask_expt & ...
        ~(ppm_axis >= water_win(1) & ppm_axis <= water_win(2));
end
% An empty artifact window is the generic-carbohydrate default: do not
% index artifact_win(1:2) or silently mask an unknown sugar's real peak.
if numel(artifact_win) >= 2
    fit_mask_expt = fit_mask_expt & ...
        ~(ppm_axis >= artifact_win(1) & ppm_axis <= artifact_win(2));
end

%% ============================================================
% Load reference curves (GISSMO simulation, real experiment) -- raw,
% not yet fit/normalized
%% ============================================================

have_gissmo = isfield(p, 'gissmo_sim_file') && ~isempty(p.gissmo_sim_file);
have_expt   = isfield(p, 'exp_file') && ~isempty(p.exp_file);

if have_gissmo
    if ~isfile(p.gissmo_sim_file)
        error('GISSMO simulated spectrum not found: %s', p.gissmo_sim_file);
    end
    raw = jsondecode(fileread(p.gissmo_sim_file));
    gissmo_ppm_raw = str2double(raw{1}); gissmo_ppm_raw = gissmo_ppm_raw(:);
    gissmo_val_raw = str2double(raw{2}); gissmo_val_raw = gissmo_val_raw(:);
    fprintf('\n===== GISSMO published simulated spectrum =====\n');
    fprintf('Loaded %s\n', p.gissmo_sim_file);
    fprintf('%d points, ppm range %.3f to %.3f\n', numel(gissmo_ppm_raw), ...
        min(gissmo_ppm_raw), max(gissmo_ppm_raw));
end

show_sucrose_diagnostics = true;
if isfield(p, 'suppress_sucrose_diagnostics') && p.suppress_sucrose_diagnostics
    show_sucrose_diagnostics = false;
end

if have_expt
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

    if ~isfield(p, 'field_label'), p.field_label = 'real'; end
    fprintf('\n===== Real %s experimental spectrum =====\n', p.field_label);
    fprintf('Read %d points from %s (expected %d)\n', numel(exp_y), p.exp_file, p.EXP_SI);

    exp_ppm_raw = p.EXP_OFFSET - (0:numel(exp_y)-1).' * ...
        (p.EXP_SW_p / p.EXP_SF) / (numel(exp_y) - 1);
    [exp_ppm_raw, ix] = sort(exp_ppm_raw);
    exp_val_raw = exp_y(ix);
    exp_val_raw = exp_val_raw - median(exp_val_raw);   % crude baseline
end

%% ============================================================
% Fit lb_Hz (+ ppm_offset) against the AUTHORITATIVE reference:
% the real experiment if available, else GISSMO's simulation, else
% just use the guess. Mirrors alanine: fit against real data first,
% theory/GISSMO is recomputed afterward at the same linewidth.
%% ============================================================

if have_expt
    fprintf('\nFitting lb_Hz + ppm_offset against the REAL experimental spectrum\n');
    if fit_receiver_phase
        fprintf('Also fitting a zero-order receiver phase for this acquisition mode.\n');
    end
    fprintf('(fit region excludes water %.2f-%.2f ppm', water_win(1), water_win(2));
    if numel(artifact_win) >= 2
        fprintf(' and the unexplained artifact band %.2f-%.2f ppm, per sucrose_refinement_notes.txt)\n', artifact_win(1), artifact_win(2));
    else
        fprintf('; no artifact band configured)\n');
    end

    if strcmp(lineshape_model, 'lorentzian')
        if fit_receiver_phase
            joint_rmse_fn = @(q) sucrose_joint_rmse_vs_ref( ...
                exp(q(1)), 0, q(2), q(3), fid_raw, dt, TD, anom_mask, fit_mask_expt, ...
                ppm_axis, exp_ppm_raw, exp_val_raw);
        else
            joint_rmse_fn = @(q) sucrose_joint_rmse_vs_ref( ...
                exp(q(1)), 0, q(2), 0, fid_raw, dt, TD, anom_mask, fit_mask_expt, ...
                ppm_axis, exp_ppm_raw, exp_val_raw);
        end
    else
        if fit_receiver_phase
            joint_rmse_fn = @(q) sucrose_joint_rmse_vs_ref( ...
                exp(q(1)), exp(q(2)), q(3), q(4), fid_raw, dt, TD, anom_mask, fit_mask_expt, ...
                ppm_axis, exp_ppm_raw, exp_val_raw);
        else
            joint_rmse_fn = @(q) sucrose_joint_rmse_vs_ref( ...
                exp(q(1)), exp(q(2)), q(3), 0, fid_raw, dt, TD, anom_mask, fit_mask_expt, ...
                ppm_axis, exp_ppm_raw, exp_val_raw);
        end
    end
    ref_label = 'experiment';
elseif have_gissmo
    fprintf('\nNo experimental file given -- fitting lb_Hz + ppm_offset against\n');
    fprintf('GISSMO''s simulation instead (fallback).\n');
    if strcmp(lineshape_model, 'lorentzian')
        if fit_receiver_phase
            joint_rmse_fn = @(q) sucrose_joint_rmse_vs_ref( ...
                exp(q(1)), 0, q(2), q(3), fid_raw, dt, TD, anom_mask, region_mask, ...
                ppm_axis, gissmo_ppm_raw, gissmo_val_raw);
        else
            joint_rmse_fn = @(q) sucrose_joint_rmse_vs_ref( ...
                exp(q(1)), 0, q(2), 0, fid_raw, dt, TD, anom_mask, region_mask, ...
                ppm_axis, gissmo_ppm_raw, gissmo_val_raw);
        end
    else
        if fit_receiver_phase
            joint_rmse_fn = @(q) sucrose_joint_rmse_vs_ref( ...
                exp(q(1)), exp(q(2)), q(3), q(4), fid_raw, dt, TD, anom_mask, region_mask, ...
                ppm_axis, gissmo_ppm_raw, gissmo_val_raw);
        else
            joint_rmse_fn = @(q) sucrose_joint_rmse_vs_ref( ...
                exp(q(1)), exp(q(2)), q(3), 0, fid_raw, dt, TD, anom_mask, region_mask, ...
                ppm_axis, gissmo_ppm_raw, gissmo_val_raw);
        end
    end
    ref_label = 'GISSMO simulation';
else
    joint_rmse_fn = [];
    ref_label = 'none';
end

if ~isempty(joint_rmse_fn)
    % Voigt fit: Lorentzian width lbL (homogeneous/T2), Gaussian width lbG
    % (inhomogeneous/shim), and the ppm calibration offset. Pure Lorentzian
    % was the old model (lbG->0); adding the Gaussian is the main lineshape
    % fix (see the Spinach<->GISSMO broadening diagnostic: 0.74->0.98 once
    % lineshape matched). lbL, lbG optimized in log-space to stay positive.
    if strcmp(lineshape_model, 'lorentzian')
        if fit_receiver_phase
            q0 = [log(p.lb_Hz_guess), 0, pi/2];
        else
            q0 = [log(p.lb_Hz_guess), 0];
        end
    else
        if fit_receiver_phase
            q0 = [log(p.lb_Hz_guess), log(0.5), 0, pi/2];
        else
            q0 = [log(p.lb_Hz_guess), log(0.5), 0];
        end
    end
    % The default fminsearch budget is too small for the Voigt + receiver
    % phase objective at some fields (900 MHz reached the default limit).
    % Increase only the optimizer budget; the objective, matrix, and bounds
    % remain unchanged so this is a convergence fix rather than a model edit.
    fit_opts = optimset('Display', 'off', 'MaxFunEvals', 5000, ...
        'MaxIter', 5000, 'TolX', 1e-8, 'TolFun', 1e-8);
    % The phase objective is periodic and can have local minima.  A single
    % zero-phase start produced visibly dispersive fits (for example +57.9
    % degrees in the sucrose overlay).  Try the standard receiver-phase
    % conventions and retain the lowest-RMSE solution.
    if fit_receiver_phase
        phase_idx = numel(q0);
        phase_seeds = [pi/2, -pi/2, 0, pi];
        best_obj = Inf;
        q_opt = q0;
        for phase_seed = phase_seeds
            q_seed = q0;
            q_seed(phase_idx) = phase_seed;
            q_try = fminsearch(joint_rmse_fn, q_seed, fit_opts);
            obj_try = joint_rmse_fn(q_try);
            if isfinite(obj_try) && obj_try < best_obj
                best_obj = obj_try;
                q_opt = q_try;
            end
        end
    else
        q_opt = fminsearch(joint_rmse_fn, q0, fit_opts);
    end
    lbL_fit = exp(q_opt(1));
    if strcmp(lineshape_model, 'lorentzian')
        lbG_fit = 0;
        ppm_offset_fit = q_opt(2);
        if fit_receiver_phase
            phase_rad_fit = atan2(sin(q_opt(3)), cos(q_opt(3)));
        else
            phase_rad_fit = 0;
        end
    else
        lbG_fit = exp(q_opt(2));
        ppm_offset_fit = q_opt(3);
        if fit_receiver_phase
            phase_rad_fit = atan2(sin(q_opt(4)), cos(q_opt(4)));
        else
            phase_rad_fit = 0;
        end
    end
    lb_Hz_fit = lbL_fit;    % kept for downstream/label continuity
    fprintf('\n===== Joint %s linewidth + calibration fit (vs %s) =====\n', upper(lineshape_model), ref_label);
    fprintf('Fitted lbL (Lorentzian) : %.4f Hz\n', lbL_fit);
    fprintf('Fitted lbG (Gaussian)   : %.4f Hz\n', lbG_fit);
    fprintf('Fitted offset           : %+.5f ppm\n', ppm_offset_fit);
    fprintf('Fitted receiver phase   : %+.2f deg\n', phase_rad_fit * 180/pi);
else
    lbL_fit = p.lb_Hz_guess;  lbG_fit = 0;
    lb_Hz_fit = lbL_fit;
    ppm_offset_fit = 0;
    phase_rad_fit = 0;
end

% Final Spinach spectrum, used everywhere downstream. Keep both normalized
% and unnormalized forms so multi-component validations can combine spectra
% before applying a common population normalization.
[spec_spinach_norm, spec_unnorm] = build_sucrose_spectrum( ...
    fid_raw, dt, lbL_fit, lbG_fit, phase_rad_fit, TD, anom_mask);

%% ============================================================
% Build final normalized reference curves at the fitted linewidth,
% and the cross-validation r/RMSE
%% ============================================================

r_sp_expt = NaN; rmse_sp_expt = NaN; exp_norm = []; exp_ppm_offset = 0;
r_sp_gissmo = NaN; rmse_sp_gissmo = NaN; gissmo_norm = []; gissmo_ppm_offset = 0;
r_gissmo_expt = NaN;

if have_expt
    exp_ppm_offset = ppm_offset_fit;
    exp_on_sp = interp1(exp_ppm_raw + exp_ppm_offset, exp_val_raw, ppm_axis, 'linear', 0);
    exp_bl = exp_on_sp - median(exp_on_sp(fit_mask_expt));
    sc_exp = max(exp_bl(anom_mask));
    if sc_exp <= 0 || ~isfinite(sc_exp), sc_exp = 1; end
    exp_norm = exp_bl / sc_exp;

    fprintf('\n===== EXPERIMENTAL ANCHOR DIAGNOSTIC =====\n');
    fprintf('Raw exp_ppm_raw range   : %.4f to %.4f ppm\n', min(exp_ppm_raw), max(exp_ppm_raw));
    fprintf('Raw exp_val_raw range   : %.6g to %.6g\n', min(exp_val_raw), max(exp_val_raw));
    [anom_exp_max, anom_exp_idx] = max(exp_bl(anom_mask));
    anom_exp_ppm = ppm_axis(anom_mask);
    fprintf('Experiment anomeric-window max : %.6g at %.4f ppm (sc_exp = %.6g)\n', ...
        anom_exp_max, anom_exp_ppm(anom_exp_idx), sc_exp);
    fprintf('exp_norm range over fit region : %.4f to %.4f  (should be ~0 to ~1)\n', ...
        min(exp_norm(fit_mask_expt)), max(exp_norm(fit_mask_expt)));
    [fitreg_max, fitreg_idx] = max(exp_norm(fit_mask_expt));
    fitreg_ppm = ppm_axis(fit_mask_expt);
    fprintf('Tallest unexplained peak inside fit region : %.4f at %.4f ppm\n', ...
        fitreg_max, fitreg_ppm(fitreg_idx));

    sp_r = spec_spinach_norm(fit_mask_expt);
    ex_r = exp_norm(fit_mask_expt);
    good = isfinite(sp_r) & isfinite(ex_r);
    r_sp_expt    = corr(sp_r(good), ex_r(good));
    rmse_sp_expt = sqrt(mean((sp_r(good) - ex_r(good)).^2));
end

if have_gissmo
    if have_expt
        % lb_Hz is already fixed (fit against experiment) -- only the
        % GISSMO-specific axis offset needs its own 1-D fit.
        gissmo_offset_fn = @(off) sucrose_joint_rmse_vs_ref( ...
            lbL_fit, lbG_fit, off, phase_rad_fit, fid_raw, dt, TD, anom_mask, region_mask, ...
            ppm_axis, gissmo_ppm_raw, gissmo_val_raw);
        gissmo_ppm_offset = fminbnd(gissmo_offset_fn, -0.05, 0.05);
    else
        gissmo_ppm_offset = ppm_offset_fit;   % already fit against GISSMO above
    end

    gissmo_on_sp = interp1(gissmo_ppm_raw + gissmo_ppm_offset, gissmo_val_raw, ppm_axis, 'linear', 0);
    gissmo_bl = gissmo_on_sp - median(gissmo_on_sp(region_mask));
    sc_gissmo = max(gissmo_bl(anom_mask));
    if sc_gissmo <= 0 || ~isfinite(sc_gissmo), sc_gissmo = 1; end
    gissmo_norm = gissmo_bl / sc_gissmo;

    sp_r = spec_spinach_norm(region_mask);
    gi_r = gissmo_norm(region_mask);
    good = isfinite(sp_r) & isfinite(gi_r);
    r_sp_gissmo    = corr(sp_r(good), gi_r(good));
    rmse_sp_gissmo = sqrt(mean((sp_r(good) - gi_r(good)).^2));

    if have_expt
        gi_r2 = gissmo_norm(fit_mask_expt);
        ex_r2 = exp_norm(fit_mask_expt);
        good2 = isfinite(gi_r2) & isfinite(ex_r2);
        r_gissmo_expt = corr(gi_r2(good2), ex_r2(good2));
    end
end

if have_expt || have_gissmo
    fprintf('\n===== Cross-validation (region %.2f-%.2f ppm, fit excludes water+artifact) =====\n', ...
        sucrose_region(1), sucrose_region(2));
    fprintf('Spinach vs experiment   r = %.4f   RMSE = %.4f\n', r_sp_expt, rmse_sp_expt);
    fprintf('Spinach vs GISSMO sim   r = %.4f   RMSE = %.4f\n', r_sp_gissmo, rmse_sp_gissmo);
    fprintf('GISSMO   vs experiment  r = %.4f\n', r_gissmo_expt);
else
    fprintf('\n===== Cross-validation skipped: no experimental or GISSMO reference supplied =====\n');
end

%% ============================================================
% With/without diagnostic: how much is the single localized mismatch
% near 3.66-3.67 ppm (likely a real per-sample shift difference in the
% matrix's near-degenerate 3.667/3.670 ppm cluster, not a pipeline bug
% -- see conversation notes) dragging down the headline numbers?
%% ============================================================

if have_expt && show_sucrose_diagnostics
    if isfield(p, 'crowded_win'), crowded_win = p.crowded_win; else, crowded_win = [3.60 3.73]; end
    fit_mask_excl_crowded = fit_mask_expt & ...
        ~(ppm_axis >= crowded_win(1) & ppm_axis <= crowded_win(2));

    sp_rc = spec_spinach_norm(fit_mask_excl_crowded);
    ex_rc = exp_norm(fit_mask_excl_crowded);
    goodc = isfinite(sp_rc) & isfinite(ex_rc);
    r_sp_expt_excl    = corr(sp_rc(goodc), ex_rc(goodc));
    rmse_sp_expt_excl = sqrt(mean((sp_rc(goodc) - ex_rc(goodc)).^2));

    fprintf('\n===== With/without the %.2f-%.2f ppm crowded-region mismatch =====\n', ...
        crowded_win(1), crowded_win(2));
    fprintf('Spinach vs experiment, WITH that region    : r = %.4f   RMSE = %.4f\n', ...
        r_sp_expt, rmse_sp_expt);
    fprintf('Spinach vs experiment, WITHOUT that region  : r = %.4f   RMSE = %.4f\n', ...
        r_sp_expt_excl, rmse_sp_expt_excl);
    fprintf('\nIf r/RMSE improve substantially once that region is excluded, the rest\n');
    fprintf('of the spectrum already agrees well and the gap is localized (consistent\n');
    fprintf('with a real per-sample shift difference in that one crowded cluster,\n');
    fprintf('not a pipeline or matrix-wide problem).\n');

    % ---- Per-peak position check on a few more ISOLATED shifts -------
    % Distinguishes two explanations for the remaining scattered
    % residual (still present even excluding the 3.66-3.73 cluster):
    %   - per-peak ppm offsets here are near-zero -> positions are fine,
    %     the residual is from linewidth/intensity mismatch instead
    %     (real samples have heterogeneous T2, unlike a clean simulation)
    %   - per-peak ppm offsets are consistently nonzero -> several
    %     shifts in the matrix differ from this sample's real values,
    %     not just the one cluster already found
    if isfield(p, 'diagnostic_peaks')
        check_peaks = p.diagnostic_peaks;
    else
        check_peaks = struct('label', {'4.207 ppm cluster','4.043 ppm cluster','3.549 ppm (Glc H2)'}, ...
            'win', {[4.17 4.24], [4.00 4.08], [3.52 3.58]});
    end
    if ~isempty(check_peaks)
        fprintf('\n===== Per-peak position check (isolated shifts) =====\n');
        for kk = 1:numel(check_peaks)
            w = check_peaks(kk).win;
            m_sp = ppm_axis >= w(1) & ppm_axis <= w(2);
            m_ex = fit_mask_expt & m_sp;
            if any(m_sp) && any(m_ex)
                tmp_ppm_sp = ppm_axis(m_sp);
                [~, i_sp] = max(spec_spinach_norm(m_sp));
                tmp_ppm_ex = ppm_axis(m_ex);
                [~, i_ex] = max(exp_norm(m_ex));
                d_ppm = tmp_ppm_ex(i_ex) - tmp_ppm_sp(i_sp);
                fprintf('%-22s Spinach %.4f ppm, experiment %.4f ppm, offset %+.4f ppm (%+.2f Hz)\n', ...
                    check_peaks(kk).label, tmp_ppm_sp(i_sp), tmp_ppm_ex(i_ex), d_ppm, d_ppm * SFO1_MHz);
            end
        end
    end
end

%% ---- Plots ----
% Optimizer/scorer callers pass p.make_plots=false to skip all figure IO
% (this function is otherwise called once per field and saves 2 PNGs; in a
% search loop that is pure overhead). Default true -> existing behavior.
make_plots = ~isfield(p,'make_plots') || isempty(p.make_plots) || p.make_plots;
if make_plots
close all;
% Field-labeled output filenames so 600/900/1100 runs don't overwrite each
% other -- the multi-field overlay set is saved as three separate PNGs.
flab = 'field';
if isfield(p,'field_label') && ~isempty(p.field_label)
    flab = regexprep(p.field_label, '\s+', '');
end
overlay_png  = sprintf('%s_vs_expt_vs_gissmo_overlay_%s.png', p.plot_prefix, flab);
residual_png = sprintf('%s_vs_expt_residual_%s.png', p.plot_prefix, flab);
% Georgia palette: Hedges (experiment), Glory Glory (Spinach), and Olympic
% (published GISSMO). Arch Black is reserved for axes and text.
col_exp     = [180 189 0] / 255;      % Hedges, #B4BD00
col_spinach = [228 0 43] / 255;       % Glory Glory, #E4002B
col_gissmo  = [0 78 96] / 255;        % Olympic, #004E60

fig1 = figure('Color','w','Position',[60 60 1600 720]);
ax1 = axes('Parent', fig1);
hold(ax1, 'on');

legend_handles = [];
if ~isempty(exp_norm)
    h_ex = plot(ax1, ppm_axis, exp_norm, '-', 'Color', col_exp, 'LineWidth', 2.6, ...
        'DisplayName', sprintf('%s experiment (%+.4f ppm calibrated)', p.field_label, exp_ppm_offset));
    legend_handles = [legend_handles h_ex];
end
h_sp = plot(ax1, ppm_axis, spec_spinach_norm, '--', 'Color', col_spinach, 'LineWidth', 2.4, ...
    'DisplayName', sprintf('Spinach FFT (%s L=%.2f/G=%.2f Hz, phase=%+.1f deg, %s)', ...
    lineshape_model, lbL_fit, lbG_fit, phase_rad_fit*180/pi, method_label));
legend_handles = [legend_handles h_sp];
if ~isempty(gissmo_norm)
    h_gi = plot(ax1, ppm_axis, gissmo_norm, ':', 'Color', col_gissmo, 'LineWidth', 2.6, ...
        'DisplayName', 'GISSMO published simulation');
    legend_handles = [legend_handles h_gi];
end

style_ax(20);
xlim(ax1, sucrose_region);
% Adaptive ylim: a fixed [-0.05 1.05] would silently clip (and hide) a
% mis-normalized curve instead of revealing it -- if anything is badly
% scaled, this should make that visually obvious rather than invisible.
y_candidates = [spec_spinach_norm(region_mask); 1];
if ~isempty(exp_norm),    y_candidates = [y_candidates; exp_norm(region_mask)];    end
if ~isempty(gissmo_norm), y_candidates = [y_candidates; gissmo_norm(region_mask)]; end
y_candidates = y_candidates(isfinite(y_candidates));
ylim(ax1, [min(-0.05, min(y_candidates)*1.1), max(1.05, max(y_candidates)*1.1)]);
xlabel(ax1, '^{1}H chemical shift (ppm)', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax1, 'Normalised intensity',        'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
title(ax1, { sprintf('%s, %s, %d-spin matrix (%s)', p.sample_label, p.field_label, nspins, method_label), ...
    sprintf('LB = %.2f Hz (fitted)    r(Spinach,expt) = %.4f    r(Spinach,GISSMO) = %.4f    r(GISSMO,expt) = %.4f', ...
    lb_Hz_fit, r_sp_expt, r_sp_gissmo, r_gissmo_expt) }, ...
    'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
legend(ax1, legend_handles, 'Location', 'northwest', 'FontSize', 16, ...
    'TextColor', 'k', 'EdgeColor', 'k', 'Color', 'w', 'LineWidth', 1.2);

exportgraphics(fig1, overlay_png, 'Resolution', 300);
fprintf('\nSaved %s\n', overlay_png);

% Export the plotted curves at native numerical resolution.  The zoom-plot
% utility uses this file rather than enlarging a small raster crop, which
% preserves smooth line shapes in standalone presentation figures.
exp_curve = nan(size(ppm_axis));
if ~isempty(exp_norm), exp_curve = exp_norm; end
gissmo_curve = nan(size(ppm_axis));
if ~isempty(gissmo_norm), gissmo_curve = gissmo_norm; end
curve_csv = sprintf('%s_curves_%s.csv', p.plot_prefix, flab);
curve_table = table(ppm_axis(:), exp_curve(:), spec_spinach_norm(:), gissmo_curve(:), ...
    'VariableNames', {'ppm', 'experiment', 'candidate_spinach', 'gissmo'});
writetable(curve_table, curve_csv);
fprintf('Saved native curve data %s\n', curve_csv);

if ~isempty(exp_norm)
    fig2 = figure('Color','w','Position',[100 100 1600 480]);
    ax2 = axes('Parent', fig2);
    hold(ax2, 'on');

    residual = spec_spinach_norm - exp_norm;
    plot(ax2, ppm_axis, residual, '-', 'Color', 'k', 'LineWidth', 1.4, ...
        'DisplayName', 'Spinach - experiment');
    yline(ax2, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');

    style_ax(18);
    xlim(ax2, sucrose_region);
    ylim_max = max(0.05, max(abs(residual(fit_mask_expt))) * 1.3);
    ylim(ax2, [-ylim_max ylim_max]);
    xlabel(ax2, '^{1}H chemical shift (ppm)', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
    ylabel(ax2, 'Residual', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
    title(ax2, sprintf('%s residual (Spinach - experiment)  |  RMSE = %.5f (water+artifact excluded from fit)', ...
        p.field_label, ...
        rmse_sp_expt), 'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');
    legend(ax2, 'Location', 'northwest', 'FontSize', 16, ...
        'TextColor', 'k', 'EdgeColor', 'k', 'Color', 'w', 'LineWidth', 1.2);

    exportgraphics(fig2, residual_png, 'Resolution', 300);
    fprintf('Saved %s\n', residual_png);
end
end   % if make_plots

%% ---- Outputs ----
results.nspins       = nspins;
results.shifts_ppm   = shifts_ppm;
results.J_Hz         = J_Hz;
results.ppm_axis     = ppm_axis;
results.spec_spinach_norm = spec_spinach_norm;
% Unnormalised spectrum is useful when combining independently simulated
% components (for example, alpha/beta anomers) before applying one common
% population/intensity normalization.
results.spec_spinach_unnorm = spec_unnorm;
results.SFO1_MHz     = SFO1_MHz;
results.O1_Hz         = O1_Hz;
results.SW_Hz         = SW_Hz;
results.B0_T          = B0_T;
results.lb_Hz_fit     = lb_Hz_fit;      % = Lorentzian part (lbL), back-compat
results.lbL_Hz        = lbL_fit;        % Voigt Lorentzian width (homogeneous)
results.lbG_Hz        = lbG_fit;        % Voigt Gaussian width (inhomogeneous)
results.lineshape_model = lineshape_model;
results.receiver_phase_deg = phase_rad_fit * 180/pi;
if have_intrinsic
    results.lb_intrinsic_Hz = lb_intrinsic;   % per-spin intrinsic FWHM used
else
    results.lb_intrinsic_Hz = [];
end
% Preserve the raw FID and acquisition timestep so callers can combine
% independently simulated components with one common linewidth/receiver
% phase before normalization (for example, Mystery Sugar alpha + beta).
results.fid_raw       = fid_raw;
results.dt            = dt;
results.phase_rad_fit = phase_rad_fit;
results.r_spinach_vs_expt    = r_sp_expt;
results.rmse_spinach_vs_expt = rmse_sp_expt;
results.exp_ppm_offset_fitted = exp_ppm_offset;
results.exp_norm = exp_norm;
results.r_spinach_vs_gissmo    = r_sp_gissmo;
results.rmse_spinach_vs_gissmo = rmse_sp_gissmo;
results.gissmo_ppm_offset_fitted = gissmo_ppm_offset;
results.gissmo_norm = gissmo_norm;
results.r_gissmo_vs_expt = r_gissmo_expt;

end

%% ============================================================
% Local helper functions
%% ============================================================

% ---- RMSE of Spinach vs ANY reference curve (GISSMO's simulation or
% ---- the real experiment -- both call this the same way), jointly as
% ---- a function of lb_Hz and a constant ppm offset applied to the
% ---- reference's axis. Used by fminsearch/fminbnd above. Mirrors
% ---- simulate_alanine_spinach_fft.m's joint_rmse_vs_expt.
function r = sucrose_joint_rmse_vs_ref( ...
    lbL, lbG, ppm_offset, phase_rad, fid_raw, dt, TD, anom_mask, ref_mask, ...
    ppm_axis, ref_ppm_raw, ref_val_raw)

    ref_on_sp_trial = interp1(ref_ppm_raw + ppm_offset, ref_val_raw, ...
        ppm_axis, 'linear', 0);
    ref_bl_trial = ref_on_sp_trial - median(ref_on_sp_trial(ref_mask));
    sc_trial = max(ref_bl_trial(anom_mask));
    if sc_trial <= 0 || ~isfinite(sc_trial), sc_trial = 1; end
    ref_norm_trial = ref_bl_trial / sc_trial;

    spec_norm = build_sucrose_spectrum(fid_raw, dt, lbL, lbG, phase_rad, TD, anom_mask);

    sp_region = spec_norm(ref_mask);
    gi_region = ref_norm_trial(ref_mask);
    good = isfinite(sp_region) & isfinite(gi_region);

    r = sqrt(mean((sp_region(good) - gi_region(good)).^2));
end

function [spec_norm, spec_unnorm] = build_sucrose_spectrum(fid_raw, dt, lbL, lbG, phase_rad, TD, anom_mask)
    t    = (0:numel(fid_raw)-1).' * dt;
    % VOIGT lineshape: a Lorentzian of FWHM lbL convolved with a Gaussian of
    % FWHM lbG. In the time domain a convolution becomes a product, so a Voigt
    % line is just an exponential decay (Lorentzian) times a Gaussian decay:
    %   Lorentzian: exp(-pi*lbL*t)
    %   Gaussian:   exp(-(pi*lbG*t)^2 / (4 ln2))
    % lbG=0 recovers the old pure-Lorentzian model. Real NMR lines are Voigt
    % (homogeneous T2 -> Lorentzian; field inhomogeneity/shim -> Gaussian), so
    % this is what lets Spinach match the experimental lineshape rather than
    % sitting at the pure-Lorentzian r~0.77.
    apod = exp(-pi * lbL * t) .* exp(-(pi * lbG * t).^2 / (4*log(2)));
    fid_apod = fid_raw .* apod;
    fid_apod(1) = fid_apod(1) / 2;

    spec_raw = fftshift(fft(fid_apod, TD));
    spec = real(spec_raw * exp(1i * phase_rad));
    spec = spec(:);
    % NOTE: do NOT fliplr() here -- see simulate_alanine_spinach_fft.m's
    % build_spinach_spectrum for why (column fliplr is a no-op, and the
    % pairing with the flipped-as-a-row ppm_axis nets out correctly).

    bl = median(spec);
    spec_bl = spec - bl;

    % Anchor to the glucose anomeric H1 peak, not global max -- a prior
    % attempt at this system normalized to global max and got burned by
    % an unexplained artifact peak being the tallest feature in real
    % data; anchoring to one specific, isolated, well-understood peak
    % avoids that failure mode even though it doesn't show up with this
    % matrix-only simulation (no artifact peak exists here).
    sc = max(spec_bl(anom_mask));
    if sc <= 0 || ~isfinite(sc), sc = 1; end

    spec_norm   = spec_bl / sc;
    spec_unnorm = spec;
end

% ---- Shared axis-formatting helper (same style as the alanine scripts) ----
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
