% Test of the implementation of the 3D FDM which can be used for solving
% PDE
% C. Weng
% DLR, Berlin
% 1st version: 14-Jun-2017
clearvars
clc


% the function to be tested and the analytic solution to its derivative
fun = @(x,y,z) (x+2).^3.*cos(pi*y).*exp(1.5*z);
dfundx = @(x,y,z) 3*(x+2).^2.*cos(pi*y).*exp(1.5*z);
dfundy = @(x,y,z) (x+2).^3.*-1*pi.*sin(pi*y).*exp(1.5*z);
dfundz = @(x,y,z) (x+2).^3.*cos(pi*y).*1.5.*exp(1.5*z);

% parameters
npx = 63;
npy = 94;
npz = 85;

xVec = linspace(-1,1,npx);
yVec = linspace(-1,1,npy);
zVec = linspace(-1,1,npz);
dx = diff(xVec([1 2]));
dy = diff(yVec([1 2]));
dz = diff(zVec([1 2]));
n = 1;  % derivative order
ooa = 8; % order of accuracy of the FDM
% mind the storage:
tic
[Dy, Dx, Dz] = getNonCompactFDmatrix3D(npy,npx,npz,dy,dx,dz,n,ooa);
toc
% Matrix free
DxF = @(f)getDerivativeMF(f,2,npy,npx,npz,dx,n,ooa);
DyF = @(f)getDerivativeMF(f,1,npy,npx,npz,dy,n,ooa);
DzF = @(f)getDerivativeMF(f,3,npy,npx,npz,dz,n,ooa);
% generate function vector
[XX,YY,ZZ] = meshgrid(xVec,yVec,zVec);
x = XX(:);
y = YY(:);
z = ZZ(:);
funVec = fun(x,y,z);
dfundxAna = dfundx(x,y,z);
dfundyAna = dfundy(x,y,z);
dfundzAna = dfundz(x,y,z);
% apply the diff. matrix
tic
dfundxNum = Dx*funVec;
dfundyNum = Dy*funVec;
dfundzNum = Dz*funVec;
toc

tic
dfundxNumMF = DxF(funVec);
dfundyNumMF = DyF(funVec);
dfundzNumMF = DzF(funVec);
toc

% error
dfundxErr = abs(dfundxNum-dfundxAna);
dfundyErr = abs(dfundyNum-dfundyAna);
dfundzErr = abs(dfundzNum-dfundzAna);
% error
dfundxErr = abs(dfundxNumMF-dfundxAna);
dfundyErr = abs(dfundyNumMF-dfundyAna);
dfundzErr = abs(dfundzNumMF-dfundzAna);
% error
dfundxErr = abs(dfundxNumMF-dfundxNum);
dfundyErr = abs(dfundyNumMF-dfundyNum);
dfundzErr = abs(dfundzNumMF-dfundzNum);
%% plot
%*****  plot the function

% dF/dx
figure(1)
subplot(311)
plot(dfundxErr)
xlabel('index'),ylabel('Error(dF/dx)')

% dF/dy
subplot(312)
plot(dfundyErr)
xlabel('index'),ylabel('Error(dF/dy)')

% dF/dz
subplot(313)
plot(dfundzErr)
xlabel('index'),ylabel('Error(dF/dz)')



%% slice
funToPlot = funVec;
figure(4)
v = x.*exp(-x.^2-y.^2-z.^2);
xslice = [-1,.0,1]; 
yslice = [1]; 
zslice = [-1 ];
sl = slice(XX,YY,ZZ,reshape(funToPlot,npy,npx,npz),xslice,yslice,zslice);
set(sl,'edgecolor',.5*[1 1 1])
xlabel('x'),ylabel('y'),zlabel('z')
colorbar







