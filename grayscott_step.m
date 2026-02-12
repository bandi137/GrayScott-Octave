function [U, V, dU_dF, dV_dF, dU_dk, dV_dk, sigmaU_map, sigmaV_map] = ...
  grayscott_step(U, V, dU_dF, dV_dF, dU_dk, dV_dk, Du, Dv, F, k, dt, ...
                 sigma_F, sigma_k, sigma_num, max_dt, clamp_max)
  % Numerically stable Gray-Scott update with uncertainty propagation.
  %
  % Stability strategy:
  % 1) split dt into explicit Euler substeps (dt_sub <= max_dt)
  % 2) clip U and V into [0, clamp_max] after every substep to prevent
  %    non-physical negative concentrations and runaway overshoots.

  n_sub = grayscott_substep_count(dt, max_dt);
  dt_sub = dt / n_sub;

  for s = 1:n_sub
    LU = grayscott_laplacian(U);
    LV = grayscott_laplacian(V);

    RU = Du*LU - U.*V.^2 + F*(1 - U);
    RV = Dv*LV + U.*V.^2 - (F + k)*V;

    U = U + RU * dt_sub;
    V = V + RV * dt_sub;

    common = V.^2;
    dRU_dF = Du*grayscott_laplacian(dU_dF) - common.*dU_dF ...
           - 2*U.*V.*dV_dF - F*dU_dF + (1 - U);
    dRV_dF = Dv*grayscott_laplacian(dV_dF) + common.*dU_dF ...
           + 2*U.*V.*dV_dF - (F + k)*dV_dF - V;
    dU_dF = dU_dF + dt_sub*dRU_dF;
    dV_dF = dV_dF + dt_sub*dRV_dF;

    dRU_dk = Du*grayscott_laplacian(dU_dk) - common.*dU_dk ...
           - 2*U.*V.*dV_dk - F*dU_dk;
    dRV_dk = Dv*grayscott_laplacian(dV_dk) + common.*dU_dk ...
           + 2*U.*V.*dV_dk - (F + k)*dV_dk - V;
    dU_dk = dU_dk + dt_sub*dRU_dk;
    dV_dk = dV_dk + dt_sub*dRV_dk;

    U = min(max(U, 0), clamp_max);
    V = min(max(V, 0), clamp_max);
  end

  sigmaU_map = sqrt((dU_dF*sigma_F).^2 + (dU_dk*sigma_k).^2 + sigma_num^2);
  sigmaV_map = sqrt((dV_dF*sigma_F).^2 + (dV_dk*sigma_k).^2 + sigma_num^2);
end
