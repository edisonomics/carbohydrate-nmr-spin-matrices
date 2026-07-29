% run_alanine_900MHz.m
%
% Driver for the 900 MHz alanine dataset -- adds a THIRD field to the
% alanine proof-of-method (600/800 are also Cece lab data). This one is
% Cece's own lab data: 900_MHz/Alanine/7 (noesypr1d, NS=16,
% 10 mM), acquired in the same 900 MHz session as the 900 sucrose
% (Sucrose/6). All pipeline logic lives in simulate_alanine_spinach_fft.m;
% this script only sets this dataset's acquisition/experimental parameters.
%
% The primary repository comparison is the lab-sample 600/800/900 set.
% The downloaded 700 MHz BMRB spectrum is not part of this repository path.

clear; close all; clc;

p.field_label   = '900 MHz';
p.output_suffix = '900MHz';

src_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fileparts(src_dir));
p.project_dir  = fullfile(repo_dir, 'outputs', 'alanine');
p.spinach_root = getenv('SPINACH_ROOT');
if isempty(p.spinach_root), error('Set SPINACH_ROOT to the Spinach installation.'); end
if ~exist(p.project_dir, 'dir'), mkdir(p.project_dir); end

% Acquisition settings (from cece_data/900_MHz/Alanine/7/acqus).
% acqus TD = 32768; the function derives ncomplex = TD/2 = 16384 acquired
% complex points and zero-fills to p.TD (= SI).
p.SFO1_MHz = 899.794229013;
p.O1_Hz    = 4229.013;
p.SW_Hz    = 9090.90909090909;
p.TD       = 32768;                 % = experimental SI (SI=32768; 2x zero-fill of 16384 complex)

% Experimental 1r (from pdata/1/procs)
p.exp_file     = fullfile(repo_dir, 'data', 'alanine', '900_MHz', '7', 'pdata', '1', '1r');
p.exp_legend_label = 'lab (900 MHz)';
p.EXP_SI      = 32768;
p.EXP_SF      = 899.789938949748;
p.EXP_SW_p    = 9090.90909090907;
p.EXP_OFFSET  = 9.81953388830781;
p.EXP_BYTORDP = 0;
p.EXP_DTYPP   = 0;
p.EXP_NC_proc = 13;

results = simulate_alanine_spinach_fft(p);
