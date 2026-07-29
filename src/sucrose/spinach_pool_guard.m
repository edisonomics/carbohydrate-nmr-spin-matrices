function info = spinach_pool_guard(requested_workers, force_reset)
%SPINACH_POOL_GUARD Keep MATLAB/Spinach pool state predictable.
%
%   spinach_pool_guard(N)       cleans a stale pool once per MATLAB session
%                               and enforces N workers when a pool exists.
%   spinach_pool_guard(N,true)  repeats the strong cleanup explicitly.
%
% Spinach's create() starts a pool when none exists, and reuses an existing
% pool even when the new system requests a different worker count.  This
% helper therefore performs the cleanup before create(), while leaving pool
% creation itself to Spinach.

if nargin < 1 || isempty(requested_workers)
    requested_workers = 1;
end
if nargin < 2 || isempty(force_reset)
    force_reset = false;
end

requested_workers = max(1, round(requested_workers));

marker = 'spinach_pool_guard_initialized';
initialized = isappdata(0, marker) && getappdata(0, marker);
if ~initialized || force_reset
    fprintf('Running Spinach pool super-clear...\n');
    if exist('smack', 'file') == 2
        % Spinach's stronger cleanup: pool, open handles, and GPU state.
        smack();
    else
        % Fallback for an incomplete/older Spinach path.
        pool = gcp('nocreate');
        if ~isempty(pool)
            delete(pool);
        end
        fclose('all');
    end
    % Root appdata survives Spinach's smack()/clear('all') call, so the
    % guard does not accidentally super-clear before every simulation.
    setappdata(0, marker, true);
end

% If another script left a pool with the wrong size, remove it.  Spinach's
% create() normally starts the requested pool, but MATLAB can race with a
% stale Processes-profile startup and silently create the profile default
% (often all available cores).  Starting the pool explicitly here makes the
% worker count deterministic; Spinach then reuses this already-correct pool.
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers ~= requested_workers
    fprintf('Replacing MATLAB pool (%d workers) with %d worker(s)...\n', ...
        pool.NumWorkers, requested_workers);
    delete(pool);
    pool = [];
end

if isempty(pool)
    try
        cluster = parcluster('local');
        pool = parpool(cluster, requested_workers);
        pool.IdleTimeout = inf;
        fprintf('Spinach pool guard started a %d-worker local pool.\n', requested_workers);
    catch err
        warning('spinach_pool_guard:pool_start_failed', ...
            'Could not pre-start a %d-worker pool (%s); Spinach will retry.', ...
            requested_workers, err.message);
        pool = gcp('nocreate');
    end
end

info = struct();
info.requested_workers = requested_workers;
info.pool_exists = ~isempty(pool);
if isempty(pool)
    info.active_workers = 0;
else
    info.active_workers = pool.NumWorkers;
end
end
