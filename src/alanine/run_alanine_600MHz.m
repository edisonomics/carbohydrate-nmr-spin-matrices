% run_alanine_600MHz.m
%
% Driver for the 600 MHz cece_data Alanine/5 dataset (noesypr1d, NS=16,
% 10 mM alanine + 10 mM DSS in D2O). All pipeline logic lives in
% simulate_alanine_spinach_fft.m -- this script only sets this
% dataset's acquisition/experimental parameters.

clear; close all; clc;

p.field_label   = '600 MHz';
p.output_suffix = '600MHz';

src_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(src_dir));
p.project_dir  = fullfile(repo_dir, 'outputs', 'alanine');
p.spinach_root = getenv('SPINACH_ROOT');
if isempty(p.spinach_root), error('Set SPINACH_ROOT to the Spinach installation.'); end
if ~exist(p.project_dir, 'dir'), mkdir(p.project_dir); end

exp_dir = fullfile(repo_dir, 'data', 'alanine', '600_MHz', '5');

% Acquisition settings (from Alanine/5/acqus)
p.SFO1_MHz = 599.764818881;
p.O1_Hz    = 2818.881;
p.SW_Hz    = 6000;
p.TD       = 32768;

% Experimental 1r (from Alanine/5/pdata/1/procs)
p.exp_file     = fullfile(exp_dir, 'pdata/1/1r');
p.exp_legend_label = '600 MHz';
p.EXP_SI      = 32768;
p.EXP_SF      = 599.761956058822;
p.EXP_SW_p    = 5999.99999999999;
p.EXP_OFFSET  = 9.775249;
p.EXP_BYTORDP = 0;
p.EXP_DTYPP   = 0;
p.EXP_NC_proc = -4;

results = simulate_alanine_spinach_fft(p);
