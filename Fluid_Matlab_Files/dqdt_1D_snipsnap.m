
function dqdt = dqdt_1D_snipsnap(t, q_int, diffmat, geom, par, bcType,dimType)

% Note that input column vector q is only interior values

switch dimType
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIMENSIONAL CASE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "Dim"
        % Collect parameters:
        alpha = par.alpha;
        Tc    = par.Tc;
        zp    = par.zp;

        % Collect geometry:
        z  = geom.z;
        rz = geom.rz;

        % Collect differentiation matrix:
        Dz = diffmat.Dz;

        % Apply boundary conditions
        switch bcType
            case "NbDt"
                % bottom Neumann
                q_l = -(Dz(1,2)*q_int(1) + Dz(1,3)*q_int(2))/Dz(1,1);
                % top Dirichlet
                q_r = 1;
            case "NbNt"
                % bottom Neumann
                q_l = -(Dz(1,2)*q_int(1) + Dz(1,3)*q_int(2))/Dz(1,1);
                % top Neumann
                q_r = -(Dz(end,end-1)*q_int(end) + Dz(end,end-2)*q_int(end-1))/Dz(end,end);
            case "DbNt"
                % bottom Dirichlet
                q_l = 1;
                % top Neumann
                q_r = -(Dz(end,end-1)*q_int(end) + Dz(end,end-2)*q_int(end-1))/Dz(end,end);
            case "DbDt"
                % bottom Dirichlet
                q_l = 1;
                % top Dirichlet
                q_r = 1;
            otherwise
                error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
        end

        %Extend to full q:
        q_full = [q_l; q_int; q_r];

        %Apply PDE to q_full:
        inv_g = 1./sqrt(rz.^2+1);
        dqdt  = inv_g.*(Dz*(inv_g.*alpha./q_full.^2 .*(Dz*q_full)))  + ...
            log(2)/Tc .* q_full .*(z<zp);

        %Chop to interior q:
        dqdt = dqdt(2:end-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NON DIMENSIONAL CASE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case "nonDim"
        % Collect parameters:
        gamma = par.gamma;
        ell   = par.ell;
        zp    = par.zp;

        % Collect geometry:
        z  = geom.z;
        rz = geom.rz;

        % Collect differentiation matrix:
        Dz = diffmat.Dz;

        % Apply boundary conditions
        switch bcType
            case "NbDt"
                % bottom Neumann
                q_l = -(Dz(1,2)*q_int(1) + Dz(1,3)*q_int(2))/Dz(1,1);
                % top Dirichlet
                q_r = ell;
            case "NbNt"
                % bottom Neumann
                q_l = -(Dz(1,2)*q_int(1) + Dz(1,3)*q_int(2))/Dz(1,1);
                % top Neumann
                q_r = -(Dz(end,end-1)*q_int(end) + Dz(end,end-2)*q_int(end-1))/Dz(end,end);
            case "DbNt"
                % bottom Dirichlet
                q_l = ell;
                % top Neumann
                q_r = -(Dz(end,end-1)*q_int(end) + Dz(end,end-2)*q_int(end-1))/Dz(end,end);
            case "DbDt"
                % bottom Dirichlet
                q_l = ell;
                % top Dirichlet
                q_r = ell;
            otherwise
                error('bcType must be one of: "NbDt", "NbNt", "DbNt", "DbDt".');
        end

        %Extend to full q:
        q_full = [q_l; q_int; q_r];

        %Apply PDE to q_full:
        inv_g = 1./sqrt(rz.^2+1);
        dqdt  = 1./gamma.*inv_g.*(Dz*(inv_g.*1./q_full.^2 .*(Dz*q_full)))  + ...
                q_full .*(z<zp);

        %Chop to interior q:
        dqdt = dqdt(2:end-1);

    otherwise
        error('dimType must be one of: "Dim", "nonDim".');
end

end