function results = run_carbohydrate_spinach_all(molecule, matrix_override)
% Run the generic Spinach model for every prepared field of a molecule.
% This is the student-facing workflow; the single-field function remains
% available for diagnostics and debugging.  An optional matrix_override
% keeps candidate-matrix validation separate from the configured GISSMO run.

if nargin < 1 || isempty(molecule), molecule = 'sucrose'; end
if nargin < 2, matrix_override = ''; end

common_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(common_dir));
summary_file = fullfile(repo_dir, 'outputs', molecule, 'prepared', ...
    'preparation_summary.csv');
if ~isfile(summary_file)
    error('Missing %s. Run prepare_carbohydrate_spectra.py first.', summary_file);
end

metadata = readtable(summary_file, 'TextType', 'string');
fields = unique(metadata.field_mhz, 'stable');
results = cell(numel(fields), 1);

for k = 1:numel(fields)
    field_label = sprintf('%g', fields(k));
    fprintf('\n===== Spinach field %s MHz (%d of %d) =====\n', ...
        field_label, k, numel(fields));
    results{k} = run_carbohydrate_spinach_field(molecule, field_label, matrix_override);
end

% Collect scalar summaries for the student-facing multifield report.
summary_rows = cell(numel(fields), 12);
for k = 1:numel(fields)
    s = results{k};
    summary_rows(k, :) = {s.field_mhz, s.field_label, s.acquisition_mode, ...
        s.pulse_program, s.r_spinach_vs_expt, s.rmse_spinach_vs_expt, ...
        s.lbL_Hz, s.lbG_Hz, s.receiver_phase_deg, s.ppm_offset_fitted, ...
        s.nspins, s.summary_file};
end
if isempty(matrix_override)
    summary_name = 'spinach_multifield_summary.csv';
else
    summary_name = 'spinach_multifield_summary_candidate.csv';
end
summary_path = fullfile(repo_dir, 'outputs', molecule, summary_name);
fid = fopen(summary_path, 'w');
if fid < 0, error('Unable to write Spinach summary: %s', summary_path); end
fprintf(fid, 'field_mhz,field_label,acquisition_mode,pulse_program,r_spinach_vs_expt,rmse_spinach_vs_expt,lbL_Hz,lbG_Hz,receiver_phase_deg,ppm_offset_fitted,nspins,summary_file\n');
for k = 1:size(summary_rows, 1)
    fprintf(fid, '%g,%s,%s,%s,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%d,%s\n', summary_rows{k, :});
end
fclose(fid);
fprintf('\nWrote %s\n', summary_path);
end
