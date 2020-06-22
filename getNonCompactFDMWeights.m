function [w,ooa,a_rats] = getNonCompactFDMWeights(dx,n,stencil)
% GETNONCOMPACTFDMWEIGHTS gets weights of the non-compact finite difference
% Inputs:
%     -dx, scalar, uniform grid spacing between each finite difference interval
%     -n, scalar, order of derivative. E.g., for d^2/dx^2 n = 2.
%     -stencil, vector, the stencil for the FDM. For example, the stencil
%               of central differences for n=2 are (cf. https://en.wikipedia.org/wiki/Finite_difference_coefficient)
%                     [-1,0,1]         (2nd order of accuray)
%                     [-2,-1,0,1,2]    (4th order of accuray)
%               the stencil of forward differences for n=3 are
%                     [0,1,2,3]        (1st order of accuray)
%                     [0,1,2,3,4]      (2nd order of accuray)
%
% Outputs:
%     -w, vector, FDM weight. The FDM estimation of u can be obtained by w*u(:)
%     -ooa, scalar, order of accuracy of the FDM
%     -a_rats, vector, character vector containing the rational approximations to the elements of the FDM coefficients
%
% C. Weng
% DLR, Berlin
% 1st version: 22-May-2017

% check size of stencil
stencil = checkInput(n,stencil);
s = length(stencil);

% generate the equation for the FDM coefficient a, i.e. sig*a=RHS
e = repmat((0:(s-1)).',1,s);
sig = repmat(stencil,s,1).^e./factorial(e);% Here it is a little different from the documentation. 
                            % I use a matlab trick to increase the precision by making sig double precision before
                            % solving the equation. I.e., the factorial(e) is divided here.
RHS = zeros(s,1);
RHS(n+1) = 1;  % due to the use of the trick, it is not factorial(n) here..
a = sig\RHS;  % got the coefficients.
w = a.'/dx^n; % got the weights.

if nargout>1
    % compute order of accuracy
    ooa = s-n;
    % compute the extra equation for special cases
    extraEq = stencil.^s./factorial(e(end,:)+1);
    tol = eps*1e4; % this is absolutely empirical...
    if abs(extraEq*a)<tol % it should be exactly zero, but there might be numerical errors...
        ooa = ooa+1;
    end
end

if nargout>2
    %     rational fraction approximations
    %     a_rat = rat(a);
    a_rats = strtrim(rats(a,30));
end

end


function stencil = checkInput(n,stencil)
stencil = stencil(:).';
s = length(stencil);
if s<=n
    error('The stencil size should be larger than n.')
end
% check lower and upper bounds of stencil
L = min(stencil);
U = max(stencil);
if L>0 || U<0 || (U-L)<1
    error('The stencil is invalid.')
end
end
