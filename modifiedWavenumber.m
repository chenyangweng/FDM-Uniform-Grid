% check the wavenumber accuracy of the numerical scheme in a numerical
% manner
%
%
% C. Weng
% DLR, Berlin
% 1st version: 22-May-2017

clear 
clc
stencil = -6:6;  % you can specify the stencil.
lambda = 1;  % wavelength
kExa = 2*pi/lambda; % exact wavenumber

% dx vector, k*dx ranges from 0 to pi
dxVec = linspace(eps,pi/kExa,50);

% get coefficients 
[a,ooa] = getNonCompactFDMWeights(1,1,stencil);
kdx = zeros(size(dxVec)); % approximated k*dx
for ind = 1:length(dxVec)
    dx = dxVec(ind); 
    x = dx*stencil;
    kdx(ind) = 1i*a*exp(-1i*kExa*x(:));
end

%% plot k*dx
figure(1)
clf
pl = plot(kExa*dxVec,kExa*dxVec,'k-',kExa*dxVec,kdx,'--');
set(pl(2),'linewidth',2)
grid minor
switch ooa
    case 1
        ooas = '1st';
    case 2
        ooas = '2nd';
    case 3
        ooas = '3rd';
    otherwise
        ooas = [num2str(ooa) 'th'];
end
legend('Exact',[ooas ' order, stencil: [' num2str(stencil) ']'])
xlabel('$k\Delta x$','interpreter','latex')
ylabel('$\tilde{k}\Delta x$','interpreter','latex')
xlim([0 pi]),ylim([0 pi])

set(findobj('type','axes'),'linewidth',0.5,'ticklength',[0.015,0.037],...
    'xminortick','off','yminortick','off','zminortick','on','xcolor','k','ycolor','k')



