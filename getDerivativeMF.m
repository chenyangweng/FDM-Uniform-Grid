function df = getDerivativeMF(f,dim,npx1,npx2,npx3,dx,n,ooa,faceIndex)
% GETDERIVATIVEMF get Matrix Free (MF) derivative of f in the dimension
% dim. MF is useful when the problem requires huge matrix storage that may
% terribly slow down the program.
% MF may be used together with MF iterative solvers, e.g. gmres
% Inputs:
%     f,   a vector, or matrix, on which the derivative is performed
%     dim, the dimension (1,2 or 3) in which the derivative is performed
%     npx1,npx2,npx3: number of points in x1, x2, x3 directions
%         respectively. Note that the condition length(f(:)) == npx1*npx2*npx3
%         must be satisfied
%     dx, spatial steps in the dim direction.
%     n, order of derivative. E.g., for d^2/dx^2 n = 2.
%     ooa, specify the order of accuracy of the FD scheme. ooa should be an
%              even number to have central difference.
%     faceIndex, if faceindex is given, then the derivative is only applied
%               to the grid specified by faceIndex. This might be helpful
%               for apply boundary conditions. Note, faceIndex could be
%               different from i_RIGHT and i_LEFT
% Output:
%     df, derivaive of f in a vector form
%
% C. Weng
% DLR, Berlin
% 1st version: 21-Jun-2017

np = npx1*npx2*npx3;
if length(f(:)) ~= np
    error('length(f(:)) must be npx1*npx2*npx3.')
end
if ~any([1 2 3]==dim)
    error('dim must be 1, 2 or 3')
end
% check ooa
if mod(ooa,2)~=0
    ooa = ooa+1;
    warning(['To have central difference, the order of accuracy is increased to ' num2str(ooa)])
end
if ~exist('faceIndex','var')
    faceIndex = [];
end
df = zeros(size(f));

% get FD coefs.
cs = getNonCompactFD_weightsForMF(dx,n,ooa);

% find the "LEFT" and "RIGHT" faces (or lines, points, depending on the
%  dimension of the problem, i.e. 3D, 2D, or 1D) in f
switch dim
    % the following i_LEFT and i_RIGHT are indices at which the FD is evaluated
    case 1
        % d^n/dx1^n
        i_LEFT = (1:npx1:np).';
        i_RIGHT = i_LEFT+(npx1-1);
        inv = 1;  % interval between the "neighboring points" in FD
    case 2
        % d^n/dx2^n
        i_LEFT = zeros(npx1*npx3,1);
        for ii = 1:npx3
            i_LEFT((ii-1)*npx1+(1:npx1)) = (ii-1)*npx1*npx2+(1:npx1).';
        end
        i_RIGHT = i_LEFT+(npx2-1)*npx1;
        inv = npx1;  % interval between the "neighboring points" in FD
    case 3
        % d^n/dx3^n
        i_LEFT = (1:1:(npx1*npx2)).';
        i_RIGHT = i_LEFT+(np-npx1*npx2);
        inv = npx1*npx2;  % interval between the "neighboring points" in FD
end
if ~isempty(faceIndex)
    i_LEFT = i_LEFT(ismember(i_LEFT,faceIndex));
    i_RIGHT = i_RIGHT(ismember(i_RIGHT,faceIndex));
end
% initialized non central index
i_NC = ones(size(df(:)));
% FD starts :)
stencil0 =0:length(cs.LEFT{1})-1; % initial stencil
% left first:
if ~isempty(i_LEFT)
    for dis = 0:length(cs.LEFT)-1  % distance from the face
        stencil = stencil0-dis;
        Fc = cs.LEFT{dis+1};   % FD coefficient
        for sind = 1:length(stencil)        % index in stencil
            s = stencil(sind);   % -1, 1 or 0 or 1 or 2...
            df(i_LEFT+inv*dis) = df(i_LEFT+inv*dis)+Fc(sind)*f(i_LEFT+inv*dis+inv*s);
        end
        %     meanwhile, update the non-central index
        i_NC(i_LEFT+inv*dis) = 0;
    end
end
% right
stencil0 = stencil0-max(stencil0);
if ~isempty(i_RIGHT)
    for dis = 0:length(cs.RIGHT)-1  % distance from the face
        stencil = stencil0+dis;
        Fc = cs.RIGHT{dis+1};   % FD coefficient
        for sind = 1:length(stencil)        % index in stencil
            s = stencil(sind);   % -1, 1 or 0 or 1 or 2...
            df(i_RIGHT-inv*dis) = df(i_RIGHT-inv*dis)+Fc(sind)*f(i_RIGHT-inv*dis+inv*s);
        end
        %     meanwhile, compute the non-central index
        i_NC(i_RIGHT-inv*dis) = 0;
    end
end
if ~isempty(faceIndex)
    i_NC = i_NC(ismember(i_NC,faceIndex));
end
% finally, the central
if ~isempty(i_NC)
    stencil = (0:(length(cs.CENTRAL)-1))-(length(cs.CENTRAL)-1)/2;
    Fc = cs.CENTRAL;   % FD coefficient
    evalInd = find(i_NC);
    for sind = 1:length(stencil)        % index in stencil
        s = stencil(sind);   % -1, 1 or 0 or 1 or 2...
        df(evalInd) = df(evalInd)+Fc(sind)*f(evalInd+inv*s);
    end
end
end

function Ceof_Struct = getNonCompactFD_weightsForMF(dx,n,ooa)
%GETNONCOMPACTFD_WEIGHTSFORMF get noncompact FD weights for Matrix Free
% (MF) computations
%   Inputs:
%     dx, spatial steps in the x direction.
%     n, order of derivative. E.g., for d^2/dx^2 n = 2.
%     ooa, specify the order of accuracy of the FD scheme. ooa should be an
%              even number to have central difference.
%   Output:
%     Ceof_Struct, coefficients of the FDM as strings.
%              Ceof_Struct.CENTRAL contains the coefficients of the central difference
%              Ceof_Struct.LEFT is a cell containing the coefficients of
%                      the left differences in coef.NC{1}, coef.NC{2}, coef.NC{3}...
%              Ceof_Struct.RIGHT is a cell containing the coefficients of
%                      the left differences in coef.NC{1}, coef.NC{2}, coef.NC{3}...
%
% C. Weng
% DLR, Berlin
% 1st version: 21-Jun-2017

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

% start with central difference
U = (sC-1)/2;
stencilC = -U:U;
[Ceof_Struct.CENTRAL,ooaCal] = getNonCompactFDMWeights(dx,n,stencilC);
if ooaCal~=ooa
    error('Order of accuracy mismatch.')
end

% the boundary treatments
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
    [w,ooaCal] = getNonCompactFDMWeights(dx,n,stencil0-ind);
    if ooaCal~=ooa
        error('Order of accuracy mismatch.')
    end
    
    % Replace the boundary points with the new weights
    if ind < mid
        % left boundary
        Ceof_Struct.LEFT{ind+1} = w;
    elseif ind > mid
        e = ceil(s-1-mid)+1;
        Ceof_Struct.RIGHT{e-ceil(ind-mid)} = w;
    end
    
end

end