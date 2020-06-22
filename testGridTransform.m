% Test of the grid transformation 
% C. Weng
% At my apartment, Lietzenburger Str., Berlin
% 1st version: 25-May-2017

% parameters controlling the distribution of the grids
ay = 1.2;
az = 1.2;
% parameters
npy = 223;
npz = 213;
np = npy*npz;
ooa = 8; % order of accuracy of the FDM

% transformation function
transfFun = @(x,a) tanh(a*x)/tanh(a);
dTransfFun = @(x,a) a*sech(a*x).^2/tanh(a);
d2TransfFun = @(x,a) -2*a^2*sech(a*x).^2.*tanh(a*x)/tanh(a);

% transfFun = @(x,a) x;
% dTransfFun = @(x,a) 1;
% d2TransfFun = @(x,a) 0;

%--------------  computational coordinate -------------
yc1D = linspace(-1,1,npy);  
zc1D = linspace(-1,1,npz);
dyc = diff(yc1D([1 2]));
dzc = diff(zc1D([1 2]));
[Dyc, Dzc] = getNonCompactFDmatrix2D(npy,npz,dyc,dzc,1,ooa);
[D2yc, D2zc] = getNonCompactFDmatrix2D(npy,npz,dyc,dzc,2,ooa);
% convert to column vector
[YYc,ZZc] = meshgrid(yc1D,zc1D);
yc = YYc(:);
zc = ZZc(:);
%--------------  physical coordinate -------------
yp = transfFun(yc,ay);
zp = transfFun(zc,az);
dyp = dTransfFun(yc,ay);
d2yp = d2TransfFun(yc,ay);
dzp = dTransfFun(zc,az);
d2zp = d2TransfFun(zc,az);
%--------------   derivatives in the physical domain, via grid
%transformation. See GridTransformation.pdf: Gathering more grid points near the boundary via grid transformation
nv = 1:np;
Dy = sparse(nv,nv,dyp.^-1)*Dyc;
D2y = sparse(nv,nv,dyp.^-2)*D2yc-sparse(nv,nv,dyp.^-3.*d2yp)*Dyc;
Dz = sparse(nv,nv,dzp.^-1)*Dzc;
D2z = sparse(nv,nv,dzp.^-2)*D2zc-sparse(nv,nv,dzp.^-3.*d2zp)*Dzc;

% generate function vector
% the function to be tested and the analytic solution to its derivative
fun = @(y,z) y.^3.*cos(pi*z);
dfundy = @(y,z) 3*y.^2.*cos(pi*z);
dfundz = @(y,z) y.^3.*-1*pi.*sin(pi*z);
d2fundy = @(y,z) 6*y.*cos(pi*z);
d2fundz = @(y,z) y.^3.*-1*pi^2.*cos(pi*z);

funVec = fun(yp,zp);
dfundyAna = dfundy(yp,zp);
dfundzAna = dfundz(yp,zp);
d2fundyAna = d2fundy(yp,zp);
d2fundzAna = d2fundz(yp,zp);

% apply the diff. matrix
dfundyNum = Dy*funVec;
dfundzNum = Dz*funVec;
d2fundyNum = D2y*funVec;
d2fundzNum = D2z*funVec;

% error
dfundyErr = abs(dfundyNum-dfundyAna);
dfundzErr = abs(dfundzNum-dfundzAna);
d2fundyErr = abs(d2fundyNum-d2fundyAna);
d2fundzErr = abs(d2fundzNum-d2fundzAna);

%% plot
YY = transfFun(YYc,ay);
ZZ = transfFun(ZZc,az);
%*****  plot the function
figure(1)
clf
surf(YY,ZZ,reshape(funVec,npz,npy))
xlabel('y'),ylabel('z'),zlabel(func2str(fun))

% dF/dy
figure(2)
subplot(221)
surf(YY,ZZ,reshape(dfundyErr,npz,npy))
shading interp 
xlabel('y'),ylabel('z'),zlabel('Error(dF/dy)')

% dF/dz
subplot(222)
surf(YY,ZZ,reshape(dfundzErr,npz,npy))
shading interp 
xlabel('y'),ylabel('z'),zlabel('Error(dF/dz)')

% d2F/dy2
subplot(223)
surf(YY,ZZ,reshape(d2fundyErr,npz,npy))
shading interp 
xlabel('y'),ylabel('z'),zlabel('Error(d^2F/dy^2)')

% d2F/dz2
subplot(224)
surf(YY,ZZ,reshape(d2fundzErr,npz,npy))
shading interp 
xlabel('y'),ylabel('z'),zlabel('Error(d^2F/dz^2)')


%%  plot the derivative in a row
inde = 1:(npz*npy);
% dF/dy
figure(3)
subplot(221)
plot(inde,dfundyAna,inde,dfundyNum,'-.')
xlabel('index'),ylabel('dF/dy')

% dF/dz
subplot(222)
plot(inde,dfundzAna,inde,dfundzNum,'-.')
xlabel('index'),ylabel('dF/dz')

% d2F/dy2
subplot(223)
plot(inde,d2fundyAna,inde,d2fundyNum,'-.')
xlabel('index'),ylabel('d^2F/dy^2')

% d2F/dz2
subplot(224)
plot(inde,d2fundzAna,inde,d2fundzNum,'-.')
xlabel('index'),ylabel('d^2F/dz^2')

%% plot mesh
fig7 = figure(7);
clf
subplot(121)
mesh(YYc,ZZc,zeros(size(YYc)))
xlabel('$\eta$','interpreter','latex'),ylabel('$\zeta$','interpreter','latex'),title('Computational grids')
view([0 90]),axis equal  tight
subplot(122)
mesh(YY,ZZ,zeros(size(YY)))
xlabel('$y$','interpreter','latex'),ylabel('$z$','interpreter','latex'),title('Physical grids')
view([0 90]),axis equal tight


set(findobj('type','axes','parent',fig7),'linewidth',0.5,'ticklength',[0.015,0.037],...
    'tickdir','out','xcolor','k','ycolor','k')




