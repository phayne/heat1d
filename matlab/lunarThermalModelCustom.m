function t = lunarThermalModelCustom( H, rhos, rhod, ks, kd, chi, ...
                                      loctime, depth, latitude, albedo )

% function t = lunarThermalModelCustom( H, rhos, rhod, ks, kd, chi, ...
%                                      loctime, depth, latitude, albedo )
% -------------------------------------------------------------------------
% Lunar thermal model for calculating temperatures using standard
% thermophysical properties.
% Inputs [dimensions]:
%   H = H-parameter in meters [1]
%   rhos = surface density in kg.m-3 [1]
%   rhod = deep density in kg.m-3 [1]
%   ks = surface conductivity in W.m-1.K-1 [1]
%   kd = conductivity at depth in W.m-1.K-1 [1]
%   chi = radiative conductivity parameter [unitless]
%   loctime = local time in (lunar) hours PAST NOON [1xn]
%   depth = depth in meters [mx1]
%   latitude = unsigned latitude in degrees [1]
%   albedo = albedo [1]
% Outputs [dimensions]:
%   t = temperature in Kelvin [mxn]
%
% Author: Paul O. Hayne (Jet Propulsion Laboratory)
% Date created: September, 2016
% Modified: September, 2017
% 
% File dependencies: Requires MATLAB "mex" file containing main thermal
% model, written in C. At the moment, this is called "heat1d_mex.c". There
% are some other dependencies provided in the header of that file.
%
% More information:
% 
% https://github.com/phayne/heat1d 
%
% Hayne, P. O., Band?eld, J. L., Siegler, M.A., Vasavada, A. R., Ghent, 
% R. R., Williams,J.-P., ? Paige, D. A. (2017). Global rego-lith 
% thermophysical properties of theMoon from the Diviner LunarRadiometer 
% Experiment. Journal ofGeophysical Research: Planets, 122,2371?2400. 
% https://doi.org/10.1002/2017JE005387
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%
% Help message: %
% %%%%%%%%%%%%%%%

if (nargin ~= 10)
    help(sprintf('%s',mfilename))
    return
end

% %%%%%%%%%%%%%%%%%
% Error checking: %
% %%%%%%%%%%%%%%%%%

% Check for existence of model files and executable:
model = 'heat1d_mex';
mexCfile = 'heat1d_mex.c';
modelHfile = 'heat1dfun.h';
orbHfile = 'orbitfun.h';
if (exist(model,'file') ~= 3)
    if (~exist(mexCfile,'file'))
        error('MATLAB mex file does not exist in this path.');
    elseif (~exist(modelHfile,'file'))
        error('Necessary file not found:\n %s',modelHfile)
    elseif (~exist(orbHfile,'file'))
        error('Necessary file not found:\n %s',orbHfile)
    else
        eval(sprintf('mex %s',mexCfile))
    end
end

% Check dimensions:
if (numel(latitude) ~= 1)
    error('Only one latitude can be accepted.');
end

if (size(loctime,1) ~= 1)
    loctime = loctime';
    if (size(loctime,1) ~= 1)
        error('Local time must be a 1-d array.');
    end
end
n = size(loctime,2);

if (size(depth,2) ~= 1)
    depth = depth';
    if (size(depth,2) ~= 1)
        error('Depth array must be a 1-d array.');
    end
end
m = size(depth, 1);

% Initialize output:
t = zeros(m,n);

% %%%%%%%%%%%%
% Constants: %
% %%%%%%%%%%%%

% rhos = 1100;    % surface density [kg.m-3]
% H = 0.06;       % H-parameter [m]
% ks = 8.1e-04;   % thermal conductivity of surface [W.m-1.K-1]
% kd = 3.3e-03;   % thermal conductivity at depth [W.m-1.K-1]
% chi = 2.7;      % radiative conductivity [unitless]

% %%%%%%%%%%%%%%%%%%%%
% Run thermal model %
% %%%%%%%%%%%%%%%%%%%%

% Optionally print time to execute program
%fprintf('Running model...\n')
%c = cputime();
[t0, z] = heat1d_mex(H, rhos, rhod, ks, kd, chi, latitude, albedo, loctime);
%fprintf('             ...done in: %f seconds\n',cputime()-c)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate results into depth array: %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:size(t,2)
    t(:,i) = interp1(z, t0(:,i), depth, 'spline');
end

% %%%%%
% END %
% %%%%%