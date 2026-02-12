function L = grayscott_laplacian(M)
  % 5-point periodic Laplacian on a unit grid.
  L = -4*M + circshift(M,[1,0]) + circshift(M,[-1,0]) ...
           + circshift(M,[0,1]) + circshift(M,[0,-1]);
end
