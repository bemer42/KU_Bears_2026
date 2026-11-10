function [Rmean, R, theta0, Rwmean, w] = avgTravelTimeFromVectorField(th, z, vth, vz, varargin)
% avgTravelTimeFromVectorField
% Compute average travel/renewal time to reach z = L by integrating
%   dtheta/dt = vth(theta,z),  dz/dt = vz(theta,z)
% using steady vector fields provided on a (theta,z) grid.
%
% REQUIRED INPUTS (sizes are enforced)
%   th   : Nth x 1 theta grid (recommended: increasing, includes 0; may include 2*pi)
%   z    : Nz  x 1 z grid (increasing; can be nonuniform)
%   vth  : Nth x Nz steady v^theta values on the grid
%   vz   : Nth x Nz steady v^z values on the grid
%
% NAME-VALUE OPTIONS
%   'z0'       : starting height (default = z(1))
%   'Nlaunch'  : number of evenly spaced launch angles in [0,2*pi) (default = 60)
%   'theta0'   : custom launch angles (overrides Nlaunch)
%   'tMax'     : max integration time (default = 1e6)
%   'RelTol'   : ODE rel tol (default = 1e-6)
%   'AbsTol'   : ODE abs tol (default = 1e-8)
%   'weights'  : weights per theta0 for weighted mean (default = [])
%
% OUTPUTS
%   Rmean  : unweighted mean travel time (NaNs omitted)
%   R      : travel times per theta0 (NaN if top not reached by tMax)
%   theta0 : launch angles used (column vector)
%   Rwmean : weighted mean (NaN if weights not provided)
%   w      : weights used (empty if none)

% ---------------- parse options ----------------
p = inputParser;
p.addParameter('z0', z(1), @(x) isnumeric(x) && isscalar(x));
p.addParameter('Nlaunch', 60, @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('theta0', [], @(x) isnumeric(x));
p.addParameter('tMax', 1e6, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('RelTol', 1e-6, @(x) isnumeric(x) && isscalar(x));
p.addParameter('AbsTol', 1e-8, @(x) isnumeric(x) && isscalar(x));
p.addParameter('weights', [], @(x) isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

% ---------------- enforce shapes ----------------
th = th(:);
z  = z(:);

Nth = numel(th);
Nz  = numel(z);

assert(isequal(size(vth), [Nth, Nz]) && isequal(size(vz), [Nth, Nz]), ...
    'Expected vth and vz to be size [numel(th) x numel(z)] (Nth x Nz).');

assert(all(diff(z) > 0), 'z must be strictly increasing.');
L = z(end);

% ---------------- periodic theta handling for interpolation ----------------
% If th includes both 0 and 2*pi, drop the last point to avoid duplication.
if abs(th(1)) < 1e-12 && abs(th(end) - 2*pi) < 1e-8
    thp  = th(1:end-1);
    vthp = vth(1:end-1,:);
    vzp  = vz(1:end-1,:);
else
    thp  = th;
    vthp = vth;
    vzp  = vz;
end

% Extend by one wrap so values are periodic at 2*pi
thExt  = [thp; thp(1) + 2*pi];
vthExt = [vthp; vthp(1,:)];
vzExt  = [vzp;  vzp(1,:)];

Fvth = griddedInterpolant({thExt, z}, vthExt, 'linear', 'nearest');
Fvz  = griddedInterpolant({thExt, z}, vzExt,  'linear', 'nearest');

% ---------------- launch angles ----------------
if isempty(opt.theta0)
    theta0 = linspace(0, 2*pi, opt.Nlaunch+1); 
    theta0(end) = [];
    theta0 = theta0(:);
else
    theta0 = opt.theta0(:);
end

% ---------------- integrate characteristics ----------------
odeOpts = odeset( ...
    'RelTol', opt.RelTol, ...
    'AbsTol', opt.AbsTol, ...
    'Events', @(t,y) hitTopEvent(t,y,L));

R = nan(numel(theta0), 1);

for k = 1:numel(theta0)
    y0 = [wrapTo2Pi(theta0(k)); opt.z0];

    rhs = @(t,y) charRHS(t,y,Fvth,Fvz);

    [~,~,tHit] = ode45(rhs, [0 opt.tMax], y0, odeOpts);

    if ~isempty(tHit)
        R(k) = tHit(1);
    else
        R(k) = NaN; % didn't reach z=L by tMax (or couldn't cross upward)
    end
end

Rmean = mean(R, 'omitnan');

% ---------------- optional weighted mean ----------------
w = opt.weights;
if isempty(w)
    Rwmean = NaN;
else
    w = w(:);
    assert(numel(w) == numel(R), 'weights must match number of theta0 launch points.');
    Rwmean = nansum(R .* w) / nansum(w);
end

end

% ================= local helper functions =================
function dydt = charRHS(~, y, Fvth, Fvz)
th = wrapTo2Pi(y(1));  % enforce periodic theta
zz = y(2);

dth = Fvth(th, zz);
dz  = Fvz(th, zz);

dydt = [dth; dz];
end

function [value, isterminal, direction] = hitTopEvent(~, y, L)
value = y(2) - L;   % event at z = L
isterminal = 1;     % stop integration
direction  = +1;    % only trigger when crossing upward
end

function th = wrapTo2Pi(th)
th = mod(th, 2*pi);
end