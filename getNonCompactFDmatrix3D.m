function [Dx1, Dx2, Dx3] = getNonCompactFDmatrix3D(npx1,npx2,npx3,dx1,dx2,dx3,n,ooa)
% GETNONCOMPACTFDMATRIX3D generates finite difference matrix for a 3D problem,
% in order to solve 3D partial differential equations.
%
% Assume the 2D matrix to be differentiated is F(x1,x2,x3), x1,x2 and x3 are
% in the 1st, 2nd and 3rd dimensions respectively. To perform the dirivative, 
% and hence solve the PDE with which F is associated, the 3D matrix F is 
% converted to a vector f with COLUMN MAJOR, i.e., f = F(:). Then, the
% derivative can be applied as Dx1*f, Dx2*f and Dx3*f
%
%Inputs:
%  npx1,npx2,npx3: number of points in x1, x2, x3 directions respectively
%  dx1,dx2,dx3: spatial steps in the x1,x2,x3 directions.
%  n: order of derivative. E.g., for d^2/dx^2 n = 2.
%  ooa: specify the order of accuracy of the FD scheme. ooa should be an
%      even number to have central difference.
%
%Outputs:
%  Dx1, Dx2, Dx3, : np*np*np sparse matrices, where np = npx1*npx2*npx3
%
% Chenyang Weng
% 1st version: 14-Jun-2017, DLR, Berlin
%
Dx1_1D = getNonCompactFDmatrix(npx1,dx1,n,ooa); 
Dx2_1D = getNonCompactFDmatrix(npx2,dx2,n,ooa);
Dx3_1D = getNonCompactFDmatrix(npx3,dx3,n,ooa);

Dx1 = kron(speye(npx3*npx2),Dx1_1D);  
Dx2 = kron(speye(npx3),kron(Dx2_1D,speye(npx1)));
Dx3 = kron(Dx3_1D,speye(npx1*npx2));

end