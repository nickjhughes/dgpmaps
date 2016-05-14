
function x = lr_left_divide(y, varargin)
%LR_LEFT_DIVIDE Left-divide y by a sum of low-rank matrices.
%
% x = lr_left_divide(y, D1, D2, ..., Dn, G1, G2, ..., Gn, H1, H2, ..., Hn, invR1, invR2, ..., invRn)
%
% Calculates x = (D1+D2+...+Dn + G1*R1*H1+G2*R2*H2+...+Gn*Rn*Hn)\y,
% where n can be any number greater than zero.

% Check inputs
if isempty(varargin)
    error('At least one D, G, H and R must be provided.');
elseif mod(length(varargin),4) ~= 0
    error('Matching D, G, H and R must be provided.');
end

% Calculate matrices
n = length(varargin)/4;
for j = 1:n
    Dj = varargin{j};
    Gj = varargin{j+n};
    Hj = varargin{j+2*n};
    invRj = varargin{j+3*n};
    if isvector(Dj)
        Dj = spdiags(Dj,0,length(Dj),length(Dj));
    end
    if ~exist('D','var')
        D = Dj;
        G = Gj;
        H = Hj;
        invR = invRj;
    else
        D = D + Dj;
        G = [G, Gj];
        H = [H; Hj];
        invR = blkdiag(invR, invRj);
    end
end

% Perform division
x = D\y;
x = H*x;
ib = D\G;
ib = invR + H*ib;
x = ib\x;
x = G*x;
x = D\(y-x);
