% Compare all 800 MHz alanine acquisitions with one fixed AX3/noesypr1d model.
% This diagnostic is deliberately separate from the single-acquisition fit.

clear; close all; clc;
src_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(src_dir));
spinach_root = getenv('SPINACH_ROOT');
data_root = fullfile(repo_dir,'data','alanine','800_MHz');
out_dir = fullfile(repo_dir,'outputs','alanine','800MHz_diagnostic');
acq_ids = {'1','8','11'};
if isempty(spinach_root) || ~isfolder(spinach_root)
    error('Set SPINACH_ROOT to the Spinach installation.');
end
if ~isfolder(out_dir), mkdir(out_dir); end
addpath(src_dir);
if exist('create','file') ~= 2, addpath(genpath(spinach_root)); end
if exist('create','file') ~= 2, error('Spinach create() not found.'); end

% Read the common acquisition axis from the 800/11 acqus file.
acqus_ref = fullfile(data_root,'11','acqus');
SFO1_MHz = bruker_param(acqus_ref,'SFO1');
O1_Hz = bruker_param(acqus_ref,'O1');
SW_Hz = bruker_param(acqus_ref,'SW_h');
TD = bruker_param(acqus_ref,'TD');
tmix_s = bruker_param(acqus_ref,'D8');
if ~isfinite(tmix_s), tmix_s = 0.05; end
fprintf('\n===== Alanine 800 MHz all-acquisition diagnostic =====\n');
fprintf('Data root: %s\n',data_root);
fprintf('SFO1 %.9f MHz, O1 %.3f Hz, SW %.3f Hz, TD %.0f, tmix %.4f s\n', ...
    SFO1_MHz,O1_Hz,SW_Hz,TD,tmix_s);

% Fixed alanine AX3 matrix: Halpha coupled equally to three methyl protons.
sys.magnet = SFO1_MHz/42.57747892;
sys.isotopes = {'1H','1H','1H','1H'};
inter.zeeman.scalar = {3.7680,1.4655,1.4655,1.4655};
inter.coupling.scalar = cell(4,4);
inter.coupling.scalar{1,2} = 7.234;
inter.coupling.scalar{1,3} = 7.234;
inter.coupling.scalar{1,4} = 7.234;
bas.formalism = 'sphten-liouv'; bas.approximation = 'none';
spin_system = basis(create(sys,inter),bas);
parameters.spins = {'1H'}; parameters.offset = O1_Hz; parameters.sweep = SW_Hz;
parameters.npoints = TD/2; parameters.zerofill = TD; parameters.axis_units = 'ppm';
parameters.invert_axis = 1; parameters.decouple = {};
parameters.rho0 = state(spin_system,'Lz','1H','cheap');
parameters.coil = state(spin_system,'L+','1H','cheap');
parameters.tmix = tmix_s; parameters.pulse_sign = 1; parameters.receiver_sign = 1;
parameters.phase_cycle = true; parameters.crusher_mode = 'none';

fprintf('Running the shared Spinach noesypr1d simulation once...\n');
tic; fid_noesy = liquid(spin_system,@alanine_noesypr1d_acquire,parameters,'nmr');
fid_noesy = fid_noesy(:);
fprintf('Shared FID complete in %.2f s; max |FID| %.6g\n',toc,max(abs(fid_noesy)));

methyl_win = [1.30 1.65]; alpha_win = [3.55 3.95];
[ppm_axis,~] = process_fid(fid_noesy,SFO1_MHz,O1_Hz,SW_Hz,TD,1.5);
score_mask = in_window(ppm_axis,methyl_win) | in_window(ppm_axis,alpha_win);
center_ppm = mean(ppm_axis(score_mask));

% Multiple phase/calibration starts are built per acquisition below.  The
% calibration grid is centered on an offset estimated from the observed
% methyl and Halpha anchors, rather than on a hard-coded field-specific
% number.  This is important for DSS-referenced data acquired on another
% instrument or with another processed OFFSET.
phase_starts = [pi/2 -pi/2 0 pi];
offset_grid = [-0.04 -0.02 0 0.02 0.04];
opts = optimset('Display','off','TolX',1e-7,'TolFun',1e-8, ...
    'MaxFunEvals',1200,'MaxIter',600);
summary_rows = cell(numel(acq_ids),10);
results = cell(numel(acq_ids),1);

for k = 1:numel(acq_ids)
    acq = acq_ids{k};
    exp_dir = fullfile(data_root,acq);
    exp_1r = fullfile(exp_dir,'pdata','1','1r');
    procs = fullfile(exp_dir,'pdata','1','procs');
    if ~isfile(exp_1r) || ~isfile(procs)
        warning('Skipping acquisition %s: missing 1r or procs.',acq);
        continue;
    end
    [exp_ppm,exp_y,proc] = read_bruker_1r_dynamic(exp_1r,procs);
    exp_norm = normalize_trace(exp_y,in_window(exp_ppm,methyl_win));
    offset_guess = estimate_reference_offset(exp_ppm,exp_norm,methyl_win,alpha_win, ...
        1.4655,3.7680);
    offset_starts = max(min(offset_guess + offset_grid,0.049),-0.049);
    starts = zeros(numel(offset_starts)*numel(phase_starts),3);
    row = 0;
    for oi = 1:numel(offset_starts)
        for pi0 = 1:numel(phase_starts)
            row = row + 1;
            starts(row,:) = [log(1.4),offset_starts(oi),phase_starts(pi0)];
        end
    end
    fprintf('Acquisition %s reference-offset guess: %+ .5f ppm (fit parameter shifts experimental axis)\n', ...
        acq,offset_guess);
    objective = @(q) fit_rmse(q,fid_noesy,SFO1_MHz,O1_Hz,SW_Hz,TD, ...
        ppm_axis,exp_ppm,exp_norm,score_mask,center_ppm);
    best_q = starts(1,:); best_val = Inf;
    for s = 1:size(starts,1)
        [q_try,val_try] = fminsearch(objective,starts(s,:),opts);
        if val_try < best_val, best_val = val_try; best_q = q_try; end
    end
    [rmse,details] = fit_rmse(best_q,fid_noesy,SFO1_MHz,O1_Hz,SW_Hz,TD, ...
        ppm_axis,exp_ppm,exp_norm,score_mask,center_ppm);
    lb_Hz = exp(best_q(1)); ppm_offset = best_q(2);
    phase_deg = wrap_to_180(best_q(3)*180/pi);
    r_fit = trace_corr(details.target,details.fit);
    [~,model_complex] = process_fid(fid_noesy,SFO1_MHz,O1_Hz,SW_Hz,TD,lb_Hz);
    model_real = real(exp(1i*best_q(3))*model_complex);
    exp_on_model = interp1(exp_ppm+ppm_offset,exp_norm,ppm_axis,'linear',0);
    X_full = [model_real(:),ones(numel(ppm_axis),1),ppm_axis(:)-center_ppm];
    model_fit = X_full*details.coef; residual = model_fit-exp_on_model;
    overlay = table(ppm_axis(:),exp_on_model(:),model_fit(:),residual(:), ...
        'VariableNames',{'ppm','experiment_norm_interp','noesypr1d_fit','residual_model_minus_expt'});
    writetable(overlay,fullfile(out_dir,sprintf('acquisition_%s_overlay.csv',acq)));
    summary_rows(k,:) = {str2double(acq),proc.NC_proc,proc.WDW,lb_Hz, ...
        ppm_offset,phase_deg,details.coef(1),r_fit,rmse,best_q(3)*180/pi};
    results{k} = struct('acquisition',acq,'ppm',ppm_axis,'experiment',exp_on_model, ...
        'fit',model_fit,'residual',residual,'lb_Hz',lb_Hz,'ppm_offset',ppm_offset, ...
        'phase_deg',phase_deg,'scale',details.coef(1),'r',r_fit,'rmse',rmse, ...
        'NC_proc',proc.NC_proc,'WDW',proc.WDW);
    fprintf('Acquisition %s: NC_proc=%g WDW=%g | LB=%.4f Hz phase=%+.2f deg scale=%+.4g | r=%.4f RMSE=%.5f\n', ...
        acq,proc.NC_proc,proc.WDW,lb_Hz,phase_deg,details.coef(1),r_fit,rmse);
end

valid = ~cellfun(@isempty,results);
summary = cell2table(summary_rows(valid,:), 'VariableNames', ...
    {'acquisition','NC_proc','WDW','lb_Hz','ppm_offset','phase_deg', ...
     'scale','r_noesypr1d_vs_expt','rmse_noesypr1d_vs_expt','raw_phase_start_deg'});
writetable(summary,fullfile(out_dir,'alanine_800MHz_acquisition_comparison.csv'));
save(fullfile(out_dir,'alanine_800MHz_acquisition_comparison.mat'), ...
    'summary','results','fid_noesy','ppm_axis','SFO1_MHz','O1_Hz','SW_Hz','TD');

fig = figure('Color','w','InvertHardcopy','off','Position',[70 70 1800 1100]);
tl = tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
title(tl,'Alanine 800 MHz: acquisition comparison with one fixed AX3/noesypr1d model', ...
    'FontSize',18,'FontWeight','bold','Color','k');
for k = 1:numel(results)
    if isempty(results{k}), continue; end
    r = results{k};
    nexttile; hold on;
    plot(r.ppm,r.experiment,'k-','LineWidth',1.5);
    plot(r.ppm,r.fit,'r--','LineWidth',1.6);
    xlim([1.1 4.1]); ylim([-0.15 1.15]); style_axis(gca,13);
    title(sprintf('Acquisition %s: r=%.4f, RMSE=%.5f',r.acquisition,r.r,r.rmse),'Color','k');
    xlabel('^1H chemical shift (ppm)'); ylabel('normalised intensity');
    legend('experiment','Spinach fit','Location','northwest');
    nexttile; plot(r.ppm,r.residual,'k-','LineWidth',1.1); hold on;
    yline(0,'-','Color',[.5 .5 .5]); xlim([1.1 4.1]); ylim([-0.18 .18]);
    style_axis(gca,13);
    title(sprintf('Acquisition %s residual; LB=%.2f Hz, phase=%+.1f deg', ...
        r.acquisition,r.lb_Hz,r.phase_deg),'Color','k');
    xlabel('^1H chemical shift (ppm)'); ylabel('fit - experiment');
end
exportgraphics(fig,fullfile(out_dir,'alanine_800MHz_acquisition_comparison.png'), ...
    'Resolution',250,'BackgroundColor','white');
fprintf('\nWrote diagnostic outputs to:\n  %s\n',out_dir);
fprintf('Choose an acquisition with physical linewidth (~1-3 Hz) and low RMSE for the next validation run.\n');

function [rmse,details] = fit_rmse(q,fid_raw,SFO1_MHz,O1_Hz,SW_Hz,TD, ...
        ppm_axis,exp_ppm,exp_norm,score_mask,center_ppm)
    lb_Hz = exp(q(1)); ppm_offset = q(2); phase_rad = q(3);
    if ~isfinite(lb_Hz) || lb_Hz<.01 || lb_Hz>30 || ...
       ~isfinite(ppm_offset) || abs(ppm_offset)>.05 || ~isfinite(phase_rad)
        rmse = 1e6; details = struct(); return;
    end
    [~,spec_complex] = process_fid(fid_raw,SFO1_MHz,O1_Hz,SW_Hz,TD,lb_Hz);
    model_real = real(exp(1i*phase_rad)*spec_complex);
    target_full = interp1(exp_ppm+ppm_offset,exp_norm,ppm_axis,'linear',0);
    m = model_real(score_mask); target = target_full(score_mask);
    x = ppm_axis(score_mask)-center_ppm; good = isfinite(m)&isfinite(target);
    m = m(good); target = target(good); x = x(good);
    X = [m(:),ones(numel(m),1),x(:)]; coef = X\target(:); fit = X*coef;
    rmse = sqrt(mean((fit-target(:)).^2));
    details = struct('coef',coef,'target',target(:),'fit',fit);
end

function [ppm_axis,spec_complex] = process_fid(fid,SFO1_MHz,O1_Hz,SW_Hz,TD,lb_Hz)
    fid = fid(:); t = (0:numel(fid)-1).' / SW_Hz;
    fid_apod = fid.*exp(-pi*lb_Hz*t); fid_apod(1) = fid_apod(1)/2;
    spec_complex = fftshift(fft(fid_apod,TD)); spec_complex = spec_complex(:);
    freq_Hz = linspace(-SW_Hz/2,SW_Hz/2,TD);
    ppm_axis = fliplr(O1_Hz/SFO1_MHz-freq_Hz/SFO1_MHz); ppm_axis = ppm_axis(:);
end

function [ppm_axis,y,p] = read_bruker_1r_dynamic(path_1r,procs_path)
    p.BYTORDP = bruker_param(procs_path,'BYTORDP');
    p.DTYPP = bruker_param(procs_path,'DTYPP');
    p.NC_proc = bruker_param(procs_path,'NC_proc');
    p.SF = bruker_param(procs_path,'SF'); p.SW_p = bruker_param(procs_path,'SW_p');
    p.OFFSET = bruker_param(procs_path,'OFFSET'); p.WDW = bruker_param(procs_path,'WDW');
    machinefmt = 'ieee-le'; if p.BYTORDP~=0, machinefmt='ieee-be'; end
    fid = fopen(path_1r,'r',machinefmt); if fid<0, error('Could not open %s',path_1r); end
    cleanup = onCleanup(@() fclose(fid));
    if p.DTYPP==0, y=fread(fid,inf,'int32'); else, y=fread(fid,inf,'double'); end
    y = double(y(:))*2^p.NC_proc; clear cleanup;
    ppm_axis = p.OFFSET-(0:numel(y)-1).'*(p.SW_p/p.SF)/(numel(y)-1);
    [ppm_axis,order] = sort(ppm_axis); y = y(order);
end

function value = bruker_param(path,name)
    value = NaN; txt = fileread(path);
    tok = regexp(txt,['##\$' regexptranslate('escape',name) '= *([^\r\n]+)'],'tokens','once');
    if isempty(tok), return; end
    raw = regexprep(strtrim(tok{1}),'<[^>]*>',''); value = str2double(strtrim(raw));
end
function y = normalize_trace(y,mask)
    y = y(:)-median(y); sc=max(y(mask)); if ~isfinite(sc)||sc<=0, sc=max(abs(y)); end
    if ~isfinite(sc)||sc<=0, sc=1; end; y=y/sc;
end
function offset_guess = estimate_reference_offset(ppm,y,methyl_win,alpha_win,expected_methyl,expected_alpha)
    methyl_mask = in_window(ppm,methyl_win); alpha_mask = in_window(ppm,alpha_win);
    offset_terms = [];
    if any(methyl_mask)
        idx = find(methyl_mask); [~,im] = max(y(methyl_mask)); observed = ppm(idx(im));
        if isfinite(observed), offset_terms(end+1) = expected_methyl-observed; end %#ok<AGROW>
    end
    if any(alpha_mask)
        idx = find(alpha_mask); [~,ia] = max(y(alpha_mask)); observed = ppm(idx(ia));
        if isfinite(observed), offset_terms(end+1) = expected_alpha-observed; end %#ok<AGROW>
    end
    if isempty(offset_terms), offset_guess = 0; else, offset_guess = median(offset_terms); end
end
function mask = in_window(x,w), mask=x>=w(1)&x<=w(2); end
function r = trace_corr(a,b)
    a=a(:)-mean(a); b=b(:)-mean(b); d=sqrt(sum(a.^2)*sum(b.^2));
    if d<=eps, r=NaN; else, r=(a'*b)/d; end
end
function d = wrap_to_180(d), d=mod(d+180,360)-180; end
function style_axis(ax,fs)
    set(ax,'Color','w','XColor','k','YColor','k','GridColor',[.83 .83 .83], ...
        'GridAlpha',.75,'FontSize',fs,'LineWidth',1.2,'Box','on','XDir','reverse','TickDir','out');
    grid(ax,'on');
end
