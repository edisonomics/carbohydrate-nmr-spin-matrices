function results = capture_sucrose_fid(field_label)
%CAPTURE_SUCROSE_FID Run the sucrose Spinach acquisition and save its FID.
%
%   results = capture_sucrose_fid('600')
%
% The default acquisition is the metadata-configured sucrose noesypr1d
% sequence.  The function returns the normal simulation results struct and
% additionally saves the raw complex FID, its time axis, and a CSV table.
% This captures the simulated time-domain signal before apodization and FFT.

if nargin < 1 || isempty(field_label)
    field_label = '600';
end
field_label = char(string(field_label));
if ~isempty(regexp(field_label, 'MHz$', 'once'))
    field_label = regexprep(field_label, 'MHz$', '');
end
if isempty(regexp(field_label, '^\d+$', 'once'))
    error('field_label must be a field such as ''600'' or ''1100''.');
end

src_dir = fileparts(mfilename('fullpath'));
addpath(src_dir);

% Run the shared, metadata-driven Spinach simulation without generating
% figures.  The returned struct contains results.fid_raw and results.dt.
results = run_sucrose_field(field_label, false);

fid = results.fid_raw(:);
dt = results.dt;
t = (0:numel(fid)-1).' * dt;

repo_dir = fileparts(fileparts(src_dir));
out_dir = fullfile(repo_dir, 'outputs', 'sucrose', [field_label 'MHz']);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

base = fullfile(out_dir, ['sucrose_' field_label 'MHz_fid']);
acquisition_mode = 'metadata-configured';
save([base '.mat'], 'fid', 't', 'dt', 'field_label', ...
    'acquisition_mode', 'results', '-v7.3');

fid_table = table(t, real(fid), imag(fid), abs(fid), ...
    'VariableNames', {'time_s', 'fid_real', 'fid_imag', 'fid_magnitude'});
writetable(fid_table, [base '.csv']);

% Plot the raw complex time-domain signal before apodization or FFT.
% The early-time panel makes the oscillatory FID visible; the full panel
% shows the complete acquisition and relaxation envelope.
fig = figure('Color', [1 1 1], 'InvertHardcopy', 'off', ...
    'Position', [100 100 1400 850]);
tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

fid_scale = max(abs(fid));
if ~isfinite(fid_scale) || fid_scale <= 0
    fid_scale = 1;
end
fid_plot = fid / fid_scale;

red = [0.80 0.00 0.00];
blue = [0.00 0.25 0.75];
black = [0.00 0.00 0.00];
grid_gray = [0.82 0.82 0.82];

ax_full = nexttile;
set(ax_full, 'Color', [1 1 1], 'XColor', [0 0 0], 'YColor', [0 0 0], ...
    'GridColor', grid_gray, 'GridAlpha', 0.8, 'FontSize', 15, ...
    'LineWidth', 1.4, 'Box', 'on');
plot(ax_full, t, real(fid_plot), '-', 'Color', red, 'LineWidth', 1.2, ...
    'DisplayName', 'real FID');
hold on;
plot(ax_full, t, imag(fid_plot), '-', 'Color', blue, 'LineWidth', 1.2, ...
    'DisplayName', 'imaginary FID');
plot(ax_full, t, abs(fid_plot), '-', 'Color', black, 'LineWidth', 1.1, ...
    'DisplayName', '|FID|');
grid(ax_full, 'on');
xlabel(ax_full, 'Time (s)', 'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold');
ylabel(ax_full, 'Normalised signal', 'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold');
title(ax_full, sprintf('Sucrose %s MHz simulated raw FID', field_label), ...
    'Color', 'k', 'FontSize', 20, 'FontWeight', 'bold');
legend(ax_full, 'Location', 'northeast', 'TextColor', 'k', ...
    'Color', 'w', 'EdgeColor', 'k', 'FontSize', 14);

ax_early = nexttile;
set(ax_early, 'Color', [1 1 1], 'XColor', [0 0 0], 'YColor', [0 0 0], ...
    'GridColor', grid_gray, 'GridAlpha', 0.8, 'FontSize', 15, ...
    'LineWidth', 1.4, 'Box', 'on');
early_end = min(t(end), 0.050);
early = t <= early_end;
plot(ax_early, t(early) * 1000, real(fid_plot(early)), '-', ...
    'Color', red, 'LineWidth', 1.4, ...
    'DisplayName', 'real FID');
hold on;
plot(ax_early, t(early) * 1000, imag(fid_plot(early)), '-', ...
    'Color', blue, 'LineWidth', 1.4, ...
    'DisplayName', 'imaginary FID');
plot(ax_early, t(early) * 1000, abs(fid_plot(early)), '-', ...
    'Color', black, 'LineWidth', 1.2, ...
    'DisplayName', '|FID|');
grid(ax_early, 'on');
xlabel(ax_early, 'Time (ms)', 'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold');
ylabel(ax_early, 'Normalised signal', 'Color', 'k', 'FontSize', 17, 'FontWeight', 'bold');
title(ax_early, 'Early-time FID detail (first 50 ms)', ...
    'Color', 'k', 'FontSize', 20, 'FontWeight', 'bold');
legend(ax_early, 'Location', 'northeast', 'TextColor', 'k', ...
    'Color', 'w', 'EdgeColor', 'k', 'FontSize', 14);

% Set these again immediately before export because MATLAB themes and
% export settings can otherwise replace an explicitly white axes face.
set([ax_full ax_early], 'Color', [1 1 1], ...
    'XColor', [0 0 0], 'YColor', [0 0 0]);
fig.InvertHardcopy = 'off';
fid_png = [base '_plot.png'];
exportgraphics(fig, fid_png, 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);

fprintf('\nSaved sucrose FID (%s MHz):\n', field_label);
fprintf('  %s\n', [base '.mat']);
fprintf('  %s\n', [base '.csv']);
fprintf('  %s\n', fid_png);
fprintf('  %d complex points, dwell time %.9g s\n', numel(fid), dt);
end
