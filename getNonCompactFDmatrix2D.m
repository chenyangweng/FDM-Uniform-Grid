function [Dx2D, Dy2D] = getNonCompactFDmatrix2D(npx,npy,dx,dy,n,ooa)
% GETNONCOMPACTFDMATRIX2D generate finite difference matrix for a 2D problem,
% in order to solve 2D partial differential equations.
%
% Assume the 2D matrix to be differentiated is F(x,y), x and y are matrix 
% "row" and "column" directions respectively. To perform the dirivative, 
% and hence solve the PDE with which F is associated, the 2D matrix F is 
% converted to a vector f with COLUMN MAJOR, i.e., f = F(:). Then, the
% derivative can be applied as Dx2D*f and Dy2D*f
%
%Inputs:
%  npx,npy: number of points in x and y directions respectively
%  dx,dy: spatial steps in the x,y directions.
%  n: order of derivative. E.g., for d^2/dx^2 n = 2.
%  ooa: specify the order of accuracy of the FD scheme. ooa should be an
%      even number to have central difference.
%
%Outputs:
%  Dx2D, Dy2D: (npx*npy)*(npx*npy) sparse matrices
%
% Chenyang Weng
% Last update: Sun Jun 30 2013, KTH, Stockholm
% Last update: 26-May-2017, DLR, Berlin
%
Dx1D = getNonCompactFDmatrix(npx,dx,n,ooa);
Dy1D = getNonCompactFDmatrix(npy,dy,n,ooa);

Dx2D = kron(Dx1D,eye(npy));
Dy2D = kron(eye(npx),Dy1D);

end