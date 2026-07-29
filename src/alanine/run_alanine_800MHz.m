% run_alanine_800MHz.m
%
% Driver for Cece's 800 MHz alanine/11 acquisition. This is deliberately
% separate from the BMRB 700 MHz driver: the 800 MHz result is from the lab
% sample and is the field used for this fit.

clear; close all; clc;

p.field_label   = '800 MHz';
p.output_suffix = '800MHz';

src_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(src_dir));
p.project_dir  = fullfile(repo_dir, 'outputs', 'alanine');
p.spinach_root = getenv('SPINACH_ROOT');
if isempty(p.spinach_root), error('Set SPINACH_ROOT to the Spinach installation.'); end
if ~exist(p.project_dir, 'dir'), mkdir(p.project_dir); end

exp_dir = fullfile(repo_dir, 'data', 'alanine', '800_MHz', '11');

% Acquisition settings from 800_MHz/alanine/11/acqus
p.SFO1_MHz = 799.713758637;
p.O1_Hz    = 3758.637;
p.SW_Hz    = 9615.385;
p.TD       = 32768;

% Processed 1r settings from 800_MHz/alanine/11/pdata/1/procs
p.exp_file     = fullfile(exp_dir, 'pdata/1/1r');
p.exp_legend_label = 'lab (800 MHz)';
p.EXP_SI      = 32768;
p.EXP_SF      = 799.71;
p.EXP_SW_p    = 9615.38461538462;
p.EXP_OFFSET  = 10.71179;
p.EXP_BYTORDP = 0;
p.EXP_DTYPP   = 0;
p.EXP_NC_proc = 10;

results = simulate_alanine_spinach_fft(p);
