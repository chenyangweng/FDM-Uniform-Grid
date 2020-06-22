function [Dx2D, Dy2D] = get2DFDMFrom1DFDM(Dx1D,Dy1D)
% GET2DFDMFROM1D converts 1D finite difference matrices to 2D, in order to
% solve 2D partial differential equations.
%
% Assume the 2D matrix to be differentiated is F(x,y), x and y are matrix 
% "row" and "column" directions respectively. To perform the dirivative, 
% and hence solve the PDE with which F is associated, the 2D matrix F is 
% converted to a vector f with COLUMN MAJOR, i.e., f = F(:). Then, the
% derivative can be applied as Dx2D*f and Dy2D*f
%
%Inputs:
%  Dx1D, Dy1D: npx*npx and npy*npy 1D finite difference matrices
%
%Outputs:
%  Dx2D, Dy2D: (npx*npy)*(npx*npy) sparse matrices
%
% Chenyang Weng
% Last update: Sun Jun 30 2013, KTH, Stockholm
% Last update: 24-May-2017, DLR, Berlin
%

npx = size(Dx1D,1);
npy = size(Dy1D,1);

Dx2D = kron(Dx1D,eye(npy));
Dy2D = kron(eye(npx),Dy1D);

end