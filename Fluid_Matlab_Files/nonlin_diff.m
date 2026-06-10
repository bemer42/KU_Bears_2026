%% Nonlinear diffusion function

function dUdt = nonlin_diff(t, U_vector, N_x, N_y, Dxx, Dyy, alpha)

U = reshape(U_vector, N_x, N_y);

D = alpha ./ (U.^2);


dU_dx = Dxx * U;
dU_dy = U * Dyy';


Flux_x = D .* dU_dx;
Flux_y = D .* dU_dy;


dFlux_dx = Dxx * Flux_x;
dFlux_dy = Flux_y * Dyy';


dU = dFlux_dx + dFlux_dy;

dUdt = dU(:);

end