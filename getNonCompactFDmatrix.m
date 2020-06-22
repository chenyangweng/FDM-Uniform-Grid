function [D,coef] = getNonCompactFDmatrix(npx,dx,n,ooa)
%   Inputs:
%     npx, number of points in the x direction.
%     dx, spatial steps in the x direction.
%     n, order of derivative. E.g., for d^2/dx^2 n = 2.
%     ooa, specify the order of accuracy of the FD scheme. ooa should be an
%              even number to have central difference.
%   Output:
%     D, sparse derivative matrices
%     coef, coefficients of the FDM as strings.
%              coef.C contains the coefficients of the central difference
%              coef.NC is a cell containing the coefficients of the non-central
%                      differences in coef.NC{1}, coef.NC{2}, coef.NC{3}...
%
% C. Weng
% DLR, Berlin
% 1st version: 22-May-2017

% check ooa
if mod(ooa,2)~=0
    ooa = ooa+1;
    warning(['To have central difference, the order of accuracy is increased to ' num2str(ooa)])
end

% generate stencil based on the ooa
s = ooa+n; % size of the sencil
sC = s;   % size of the sencil for central difference
if mod(n,2)==0
    % if n is an even number, the calculation of ooa for the central differnce is a special case;
    sC = sC-1;  % As you may see, Sc is always an odd number
end

% examine the size of npx
if npx<s
    error('Not enough grid points for the specified order of accuracy.')
end

% start with central difference
U = (sC-1)/2;
stencilC = -U:U;
[wC,ooaCal,a_rats] = getNonCompactFDMWeights(dx,n,stencilC);
if ooaCal~=ooa
    error('Order of accuracy mismatch.')
end
if nargout>1
    % the user is interested in the FDM coefficients...
    coef.C = a_rats;
end
% Put up the central differences
vix = ones(npx,1);

D = spdiags(kron(wC,vix),stencilC,npx,npx);

% now, the boundary treatments
if nargout>1
    % the user is interested in the FDM coefficients...
    nnc = floor(s/2)*2; % number of non-central differences
    coef.NC{nnc} = [];
    inc = 1;
end
stencil0 = 0:(s-1); % inital stencil
for ind = stencil0
    % determine whether the weight w belongs to the left or the right
    % boundary
    mid = (s-1)/2; % mid points
    if ind == mid
        % This happens when s is an odd number
        % we should skip this case because there is no need to compute the weight
        % for the central differnce.
        % If s is an even number, then (s-1)/2 is a non-integer and the
        % condition ind == mid won't be satisfied.
        continue
    end
    % get weights
    [w,ooaCal,a_rats] = getNonCompactFDMWeights(dx,n,stencil0-ind);
    if ooaCal~=ooa
        error('Order of accuracy mismatch.')
    end
    
    % Replace the boundary points with the new weights
    if ind < mid
        % left boundary
        D(ind+1,1:s) = w;
    elseif ind > mid
        D(end-(s-ind-1),(end-s+1):end) = w;
    end
    if nargout>1
        % the user is interested in the FDM coefficients...
        coef.NC{inc} = a_rats;
        inc = inc+1;
    end
end
end

