1;

%!test
%! M = ones(8, 8);
%! L = grayscott_laplacian(M);
%! assert(max(abs(L(:))), 0, 1e-12);

%!test
%! assert(grayscott_substep_count(1.0, 0.25), 4);
%! assert(grayscott_substep_count(0.1, 0.25), 1);

%!test
%! n = 20;
%! U = ones(n, n);
%! V = zeros(n, n);
%! U(8:12, 8:12) = 0.5;
%! V(8:12, 8:12) = 0.25;
%! z = zeros(n, n);
%! [U2, V2, dU_dF, dV_dF, dU_dk, dV_dk, sU, sV] = ...
%!   grayscott_step(U, V, z, z, z, z, 0.16, 0.08, 0.03, 0.057, 1.0, ...
%!                 5e-4, 5e-4, 1e-6, 0.25, 1.2);
%! assert(all(isfinite(U2(:))));
%! assert(all(isfinite(V2(:))));
%! assert(min(U2(:)) >= 0);
%! assert(min(V2(:)) >= 0);
%! assert(max(U2(:)) <= 1.2 + 1e-12);
%! assert(max(V2(:)) <= 1.2 + 1e-12);
%! assert(all(isfinite(sU(:))));
%! assert(all(isfinite(sV(:))));
%! assert(norm(dV_dk(:), 2) > 0);
