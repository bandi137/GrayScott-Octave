function n_sub = grayscott_substep_count(dt, max_dt)
  % Split a user time step into stable explicit-Euler substeps.
  n_sub = max(1, ceil(dt / max_dt));
end
