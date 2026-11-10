function Renewal_Time = ND_compute_Renewal_Time(z_0, th_0_vec, Vth, Vz_hat, Tc, geom_hat)

% -------------------------
% Collect geometry and parameters: 
% -------------------------
th    = geom_hat.th; 
z_hat = geom_hat.z_hat;
S_hat = geom_hat.S_hat;

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
    th_ext = [th; 2*pi];
    Uth    = ex(Vth) ./ ex(S_hat);
    Uz_hat = ex(Vz_hat) ./ ex(S_hat);

    % -------------------------
    % Create gridded interpolants:
    % -------------------------
    Fth = griddedInterpolant({th_ext, z_hat}, Uth,    "linear", "nearest");
    Fz  = griddedInterpolant({th_ext, z_hat}, Uz_hat, "linear", "nearest");

    % -------------------------
    % Define ODE:
    % -------------------------
    rhs = @(t,y) rhsPeriodic(t, y, Fth, Fz, thMin, P);

    % -------------------------
    % Event to stop when z reaches L
    % -------------------------    
    opts = odeset("Events", @(t,y) eventTop(t,y), ...
        "RelTol", 1e-7, "AbsTol", 1e-9);

    % -------------------------
    % Integrate:
    % -------------------------    
    tSpan = [0, log(2)/Tc*1e6];  
    y0    = [th_0_vec(i); z_0];
    [t,y,te,ye,ie] = ode45(rhs, tSpan, y0, opts);

    % -------------------------
    % System check:
    % -------------------------
    if isempty(te)
        warning("Trajectory did not reach z = 1 within tSpan (or Vz never becomes positive enough).");
        travelTime = NaN;
    else
        travelTime = te(end);
    end

    % -------------------------
    % Store renewal time:
    % -------------------------
    Time(i) = travelTime;

end

Renewal_Time = (Tc/log(2))*mean(Time)/24;

end

%% Helpers (local functions)
function [value, isterminal, direction] = eventTop(~, y)
value      = y(2) - 1;
isterminal = 1;
direction  = +1;
end

function dydt = rhsPeriodic(~, y, Fth, Fz, thMin, P)
th = thMin + mod(y(1) - thMin, P);
z  = y(2);
dydt = [Fth(th,z); Fz(th,z)];
end