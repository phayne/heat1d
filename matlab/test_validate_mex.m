% test_validate_mex.m
% Validation test suite for the heat1d MEX thermal model.
%
% Runs standard test cases and prints PASS/FAIL results matching the
% style used by the C and Python validation suites.
%
% Usage: test_validate_mex
%
% Author: Paul Hayne / auto-generated

%% Setup: compile MEX if needed
if exist('heat1d_mex','file') ~= 3
    fprintf('Compiling heat1d_mex.c ...\n');
    mex heat1d_mex.c
end

%% Standard lunar regolith parameters (Hayne et al., 2017)
H      = 0.06;       % H-parameter [m]
rhos   = 1100;        % surface density [kg/m^3]
rhod   = 1800;        % density at depth [kg/m^3]
ks     = 7.4e-4;      % surface conductivity [W/m/K]
kd     = 3.4e-3;      % conductivity at depth [W/m/K]
chi    = 2.7;         % radiative conductivity parameter
albedo = 0.12;        % surface albedo

% Local time output array (hours past noon, 0-24)
loctime = linspace(0, 24, 481);
loctime = loctime(1:end-1);  % 480 points, open interval

% Depth array for output interpolation
depth = logspace(-3, 0, 50)';

n_pass = 0;
n_fail = 0;

fprintf('\n=== heat1d MEX Validation Suite ===\n\n');

%% Test 1: Equator peak temperature
lat_deg = 0;
solver = 0;  % explicit
[t0, z] = heat1d_mex(H, rhos, rhod, ks, kd, chi, lat_deg, albedo, loctime, solver);
Tsurf = t0(1,:);
Tmax = max(Tsurf);
Tmin = min(Tsurf);
Tmean = mean(Tsurf);

% Hayne et al. (2017) Table A2: equator peak ~385 K, min ~95 K
tol_peak = 10;  % K tolerance
pass_peak = abs(Tmax - 385) < tol_peak;
if pass_peak
    n_pass = n_pass + 1;
    fprintf('  [PASS] equator_Tmax  Tmax=%.1f K (ref: 385 K, tol: %d K)\n', Tmax, tol_peak);
else
    n_fail = n_fail + 1;
    fprintf('  [FAIL] equator_Tmax  Tmax=%.1f K (ref: 385 K, tol: %d K)\n', Tmax, tol_peak);
end

tol_min = 15;  % K tolerance (nighttime harder to match exactly)
pass_min = abs(Tmin - 95) < tol_min;
if pass_min
    n_pass = n_pass + 1;
    fprintf('  [PASS] equator_Tmin  Tmin=%.1f K (ref: 95 K, tol: %d K)\n', Tmin, tol_min);
else
    n_fail = n_fail + 1;
    fprintf('  [FAIL] equator_Tmin  Tmin=%.1f K (ref: 95 K, tol: %d K)\n', Tmin, tol_min);
end

%% Test 2: lat=0 does not crash (already ran above, but explicitly note it)
n_pass = n_pass + 1;
fprintf('  [PASS] lat_zero_nocrash  lat=0 ran without error\n');

%% Test 3: Solver consistency (explicit vs CN vs implicit)
solver_names = {'explicit', 'crank-nicolson', 'implicit'};
solver_ids   = [0, 1, 2];
Tmax_all = zeros(1,3);
solver_times = zeros(1,3);

for s = 1:3
    tic;
    [ts, ~] = heat1d_mex(H, rhos, rhod, ks, kd, chi, lat_deg, albedo, loctime, solver_ids(s));
    solver_times(s) = toc;
    Tmax_all(s) = max(ts(1,:));
end

% All solvers should agree on Tmax within 1 K
tol_solver = 1.0;
Tmax_range = max(Tmax_all) - min(Tmax_all);
pass_solver = Tmax_range < tol_solver;
if pass_solver
    n_pass = n_pass + 1;
    fprintf('  [PASS] solver_consistency  Tmax range=%.2f K across solvers (tol: %.1f K)\n', ...
        Tmax_range, tol_solver);
else
    n_fail = n_fail + 1;
    fprintf('  [FAIL] solver_consistency  Tmax range=%.2f K across solvers (tol: %.1f K)\n', ...
        Tmax_range, tol_solver);
end

for s = 1:3
    fprintf('    %16s: Tmax=%.1f K\n', solver_names{s}, Tmax_all(s));
end

%% Test 4: Implicit solver stability (no NaN or negative temps)
[ti, ~] = heat1d_mex(H, rhos, rhod, ks, kd, chi, 45, albedo, loctime, 2);
pass_stable = ~any(isnan(ti(:))) && all(ti(:) > 0);
if pass_stable
    n_pass = n_pass + 1;
    fprintf('  [PASS] implicit_stability  lat=45, no NaN or negative temps\n');
else
    n_fail = n_fail + 1;
    fprintf('  [FAIL] implicit_stability  lat=45, NaN or negative temps found\n');
end

%% Test 5: Multi-latitude monotonic peak temperature
latitudes = [0, 15, 30, 45, 60, 75];
Tmax_lat = zeros(size(latitudes));
for i = 1:length(latitudes)
    [tl, ~] = heat1d_mex(H, rhos, rhod, ks, kd, chi, latitudes(i), albedo, loctime, 0);
    Tmax_lat(i) = max(tl(1,:));
end
% Peak temperature should decrease with latitude
pass_mono = all(diff(Tmax_lat) < 0);
if pass_mono
    n_pass = n_pass + 1;
    fprintf('  [PASS] latitude_monotonic  Tmax decreases with latitude\n');
else
    n_fail = n_fail + 1;
    fprintf('  [FAIL] latitude_monotonic  Tmax does not decrease monotonically\n');
end
for i = 1:length(latitudes)
    fprintf('    lat=%2d: Tmax=%.1f K\n', latitudes(i), Tmax_lat(i));
end

%% Summary
fprintf('\n--- Summary ---\n');
fprintf('  %d passed, %d failed\n\n', n_pass, n_fail);

%% Solver speed comparison
fprintf('  Solver speed comparison (equator):\n');
for s = 1:3
    fprintf('    %16s: %.3f s\n', solver_names{s}, solver_times(s));
end
fprintf('\n');
