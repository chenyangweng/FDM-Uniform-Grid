% check the convergence of the numerical scheme for an acoustic boundary
% layer problem
%
%
% C. Weng
% DLR, Berlin
% 1st version: 23-May-2017

clear 
clc

% define parameters
Sh = 100; % shear wavenumber


%-----   computing the error
NVec = round(logspace(1.5,3,15));
err(length(NVec)) = 0;  % initialize the error vector
ooa = 6;            % order of accuracy
for ind = 1:length(NVec)
    N = NVec(ind);
    y = linspace(0,1,N);
    dy = diff(y([1 2]));
    D2 = getNonCompactFDmatrix(N,dy,2,ooa);  % derivative matrix
    % construct the eq.
    LHS = D2-1i*Sh^2*speye(N);
    % boundary condition
    LHS([1 end],:) = 0;
    LHS(1,1) = 1;
    LHS(end,end) = 1;
    RHS = zeros(N,1);
    RHS(1) = 1;
    v = LHS\RHS;
    vEx = exp(-sqrt(1i)*Sh*y(:));
    err(ind) = max(abs(v-vEx));
end


%% fit the error to polynomials
polyErr = (1./NVec).^ooa;
polyErr = polyErr/polyErr(1)*2;
%------------- plot errors
figure(1)
clf
pl = loglog(NVec,err,'-',NVec,polyErr,'k--');

xlabel('$N$','interpreter','latex'),
ylabel('$\mathcal{E}$','interpreter','latex')
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
set(pl(1),'displayname',[ooas ' order'],'linewidth',1)
legend(pl(1))
set(gca,'linewidth',0.5,'ticklength',[0.025,0.037],'xtick',10.^(0:5),'ytick',10.^(-16:0),...
    'xminortick','on','yminortick','on','zminortick','on','xcolor','k','ycolor','k')

return
%% --- plot Asymptotic solution

% plot exact solution
fig2 = figure(2);
clf
subplot(121)
semilogx(y,real(vEx),'k-',y,real(v),'r.')
legend('Asymptotic',['Numeric, \itN\rm=' num2str(N)])
xlabel('$y$','interpreter','latex'),
ylabel('$\Re({v}/{v_{\mbox{w}}})$','interpreter','latex'),
subplot(122)
semilogx(y,imag(vEx),'k-',y,imag(v),'r.')
legend('Asymptotic',['Numeric, \itN\rm=' num2str(N)])
xlabel('$y$','interpreter','latex'),
ylabel('$\Im({v}/{v_{\mbox{w}}})$','interpreter','latex')
set(findobj('type','axes','parent',fig2),'linewidth',0.5,'ticklength',[0.025,0.037],'xtick',10.^(-10:3),...
    'xminortick','on','yminortick','off','zminortick','on','xcolor','k','ycolor','k')



