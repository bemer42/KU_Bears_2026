%% Nonlinear diffusion function

function dUdt = nonlin_diff(t, U_vector, N_x, N_y, Dx, Dy, alpha)

U = reshape(U_vector, N_x, N_y);

D = alpha ./ (U.^2);


dU_dx = Dx * U;
dU_dy = U * Dy';


Flux_x = D .* dU_dx;
Flux_y = D .* dU_dy;


dFlux_dx = Dx * Flux_x;
dFlux_dy = Flux_y * Dy';


dU = dFlux_dx + dFlux_dy;

dUdt = dU(:);

end