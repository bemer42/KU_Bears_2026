function Renewal_Time = compute_Renewal_Time(z_0,th_0_vec, Vth, Vz, geom)

% -------------------------
% Collect geometry: 
% -------------------------
th  = geom.th; 
z   = geom.z;
G   = geom.G;
L   = geom.L; 

% -------------------------
% Initial condition vector:
% -------------------------
N_vec    = length(th_0_vec);

% -------------------------
% Initialize renewal time:
% -------------------------
Time = zeros(N_vec,1);

% -------------------------
% Loop through intial conditions:
% -------------------------
ex = @(x) [x; x(1,:)];
for i = 1:N_vec

    % -------------------------
    % Period parameters:
    % -------------------------
    thMin = 0;
    P     = 2*pi;

    % -------------------------
    % Extend theta grid to include 2*pi and extend fields:
    % -------------------------
    th_ext  = [th; 2*pi];
    Uth = ex(Vth) ./ ex(sqrt(G));
    Uz  = ex(Vz) ./ ex(sqrt(G));

    % -------------------------
    % Create gridded interpolants:
    % -------------------------
    Fth = griddedInterpolant({th_ext, z}, Uth, "linear", "nearest");
    Fz  = griddedInterpolant({th_ext, z}, Uz,  "linear", "nearest");

    % -------------------------
    % Define ODE:
    % -------------------------
    rhs = @(t,y) rhsPeriodic(t, y, Fth, Fz, thMin, P);

    % -------------------------
    % Event to stop when z reaches L
    % -------------------------
    opts = odeset("Events", @(t,y) eventTop(t,y,L), ...
        "RelTol", 1e-7, "AbsTol", 1e-9);

    % -------------------------
    % Integrate:
    % -------------------------
    tSpan = [0, 1e3];
    y0    = [th_0_vec(i); z_0];
    [t,y,te,ye,ie] = ode45(rhs, tSpan, y0, opts);

    % -------------------------
    % System check:
    % -------------------------
    if isempty(te)
        warning("Trajectory did not reach z = L within tSpan (or Vz never becomes positive enough).");
        travelTime = NaN;
    else
        travelTime = te(end);
    end

    % -------------------------
    % Store renewal time:
    % -------------------------
    Time(i) = travelTime;

end

Renewal_Time = mean(Time)/24;

end

%% Helpers (local functions)
function [value, isterminal, direction] = eventTop(~, y, L)
value      = y(2) - L;
isterminal = 1;
direction  = +1;
end

function dydt = rhsPeriodic(~, y, Fth, Fz, thMin, P)
th = thMin + mod(y(1) - thMin, P);
z  = y(2);
dydt = [Fth(th,z); Fz(th,z)];
end