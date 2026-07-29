% Optional 800 MHz condition run. This dataset is kept separate from the
% primary 600/900/1100 matrix-validation set because its sample metadata differ.
clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));
run_sucrose_field('800');
