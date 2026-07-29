function fid = alanine_noesypr1d_acquire(spin_system, parameters, H, R, K)
%ALANINE_NOESYPR1D_ACQUIRE Bruker noesypr1d preparation plus acquisition.
%
% This follows the pulse program used by cece_data/600_MHz/Alanine/5:
%
%   d1 presat, p1 ph1, p1 ph2, d8 presat/mixing, p1 ph3, acquire ph31
%
% Presaturation is optional. When parameters.presat_nu1_Hz is supplied,
% the d1 and/or tmix delays evolve under a weak continuous RF field at the
% transmitter offset. This lets us test the real noesypr1d water-presat
% physics against alanine without inventing per-peak response factors.

grumble(spin_system, parameters, H, R, K);

L = H + 1i*R + 1i*K;

if isfield(parameters, 'decouple')
    [L, parameters.rho0] = decouple(spin_system, L, parameters.rho0, parameters.decouple);
end

Lx = operator(spin_system, 'Lx', parameters.spins{1});
Ly = operator(spin_system, 'Ly', parameters.spins{1});

pulse_sign = get_scalar(parameters, 'pulse_sign', 1);
receiver_sign = get_scalar(parameters, 'receiver_sign', 1);
final_pulse_angle = get_scalar(parameters, 'final_pulse_angle', pi/2);
tmix = get_scalar(parameters, 'tmix', 0.05);
d1 = get_scalar(parameters, 'd1', 0);
presat_nu1_Hz = get_scalar(parameters, 'presat_nu1_Hz', 0);

presat_d1 = false;
if isfield(parameters, 'presat_d1') && ~isempty(parameters.presat_d1)
    presat_d1 = logical(parameters.presat_d1);
end

presat_tmix = false;
if isfield(parameters, 'presat_tmix') && ~isempty(parameters.presat_tmix)
    presat_tmix = logical(parameters.presat_tmix);
end

phase_cycle = true;
if isfield(parameters, 'phase_cycle') && ~isempty(parameters.phase_cycle)
    phase_cycle = logical(parameters.phase_cycle);
end

crusher_mode = 'none';
if isfield(parameters, 'crusher_mode') && ~isempty(parameters.crusher_mode)
    crusher_mode = lower(strtrim(parameters.crusher_mode));
end

if isfield(parameters, 'coil') && ~isempty(parameters.coil)
    coil = parameters.coil;
else
    coil = state(spin_system, 'L+', parameters.spins{1}, 'cheap');
end

if presat_nu1_Hz > 0
    L_presat = L + 2*pi*presat_nu1_Hz*Lx;
else
    L_presat = L;
end

% Bruker phase codes: 0=x, 1=y, 2=-x, 3=-y.
if phase_cycle
    ph1  = [0 2];
    ph2  = [0 0 0 0 0 0 0 0 2 2 2 2 2 2 2 2];
    ph3  = [0 0 2 2 1 1 3 3];
    ph31 = [0 2 2 0 1 3 3 1 2 0 0 2 3 1 1 3];
else
    ph1 = 0; ph2 = 0; ph3 = 0; ph31 = 0;
end

rho_stack = repmat(parameters.rho0, 1, numel(ph31));

for n = 1:numel(ph31)
    rho = parameters.rho0;

    if d1 > 0
        if presat_d1 && presat_nu1_Hz > 0
            rho = evolution(spin_system, L_presat, [], rho, d1, 1, 'final');
        else
            rho = evolution(spin_system, L, [], rho, d1, 1, 'final');
        end
    end

    [op, ang] = phase_pulse(Lx, Ly, ph1(wrap_index(n, numel(ph1))), pulse_sign*pi/2);
    rho = step(spin_system, op, rho, ang);

    [op, ang] = phase_pulse(Lx, Ly, ph2(wrap_index(n, numel(ph2))), pulse_sign*pi/2);
    rho = step(spin_system, op, rho, ang);

    rho_stack(:, n) = rho;
end

if tmix > 0
    if presat_tmix && presat_nu1_Hz > 0
        rho_stack = evolution(spin_system, L_presat, [], rho_stack, tmix, 1, 'final');
    else
        rho_stack = evolution(spin_system, L, [], rho_stack, tmix, 1, 'final');
    end
end

switch crusher_mode
    case 'none'
    case 'zeroq'
        rho_stack = coherence(spin_system, rho_stack, {{parameters.spins{1}, 0}});
    otherwise
        error('Unknown crusher_mode "%s"; use "none" or "zeroq".', crusher_mode);
end

for n = 1:numel(ph31)
    [op, ang] = phase_pulse(Lx, Ly, ph3(wrap_index(n, numel(ph3))), pulse_sign*final_pulse_angle);
    rho_stack(:, n) = step(spin_system, op, rho_stack(:, n), ang);
end

fid_stack = evolution(spin_system, L, coil, rho_stack, ...
    1/parameters.sweep, parameters.npoints-1, 'observable');

weights = exp(1i*receiver_sign*pi*ph31(:)/2);
if size(fid_stack, 2) == numel(ph31)
    fid = fid_stack * weights;
elseif size(fid_stack, 1) == numel(ph31)
    fid = (weights.' * fid_stack).';
else
    error('Unexpected FID stack size %s for %d phase-cycle steps.', ...
        mat2str(size(fid_stack)), numel(ph31));
end

fid = fid / numel(ph31);
fid = fid(:);

end

function idx = wrap_index(n, len)
idx = mod(n-1, len) + 1;
end

function value = get_scalar(parameters, field_name, default_value)
value = default_value;
if isfield(parameters, field_name) && ~isempty(parameters.(field_name))
    value = parameters.(field_name);
end
end

function [op, angle] = phase_pulse(Lx, Ly, phase_code, base_angle)
switch phase_code
    case 0
        op = Lx; angle = base_angle;
    case 1
        op = Ly; angle = base_angle;
    case 2
        op = Lx; angle = -base_angle;
    case 3
        op = Ly; angle = -base_angle;
    otherwise
        error('Unknown phase code %d.', phase_code);
end
end

function grumble(spin_system, parameters, H, R, K)
if ~ismember(spin_system.bas.formalism, {'sphten-liouv', 'zeeman-liouv'})
    error('alanine_noesypr1d_acquire requires a Liouville-space formalism.');
end
if (~isnumeric(H)) || (~isnumeric(R)) || (~isnumeric(K)) || ...
   (~ismatrix(H)) || (~ismatrix(R)) || (~ismatrix(K))
    error('H, R and K must be matrices.');
end
if (~all(size(H) == size(R))) || (~all(size(R) == size(K)))
    error('H, R and K must have the same dimensions.');
end
if ~isfield(parameters, 'sweep') || numel(parameters.sweep) ~= 1 || parameters.sweep <= 0
    error('parameters.sweep must be a positive scalar.');
end
if ~isfield(parameters, 'npoints') || numel(parameters.npoints) ~= 1 || parameters.npoints < 1
    error('parameters.npoints must be a positive scalar.');
end
if ~isfield(parameters, 'spins') || numel(parameters.spins) ~= 1
    error('parameters.spins must contain one isotope label.');
end
if ~isfield(parameters, 'rho0')
    error('parameters.rho0 is required.');
end
end
