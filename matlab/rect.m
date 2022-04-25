function [xx,yy,zz] = rect(L, H, P, U, W)

% RECT  Generate rectangle.
%
%   [X,Y,Z] = RECT(L,H,P,U,W) generates three 2-by-2 matrices so
%   that SURF(X,Y,Z) produces a rectangle of length L and height
%   H, centered at the point P, and laid in the plan defined by
%   the vectors U and W in such a way that its length stretches
%   out in the direction of U.
%
%   [X,Y,Z] = RECT uses L = 1, H = 1, P = [0 0 0], U = [1 0 0],
%   and W = [0 1 0].
%
%   The following forms are permitted (note: missing arguments
%   are set as above):
%   [X,Y,Z] = RECT(L,H,P)
%   [X,Y,Z] = RECT(L,H)
%   [X,Y,Z] = RECT(L)
%
%   RECT and RECT(...) graph the rectangle as a SURFACE
%   and do not return anything.
%
%   See also SPHERE, DISK, CYLINDER.

%   Charlles R. A. Abreu, Escola de Quimica, Universidade
%   Federal do Rio de Janeiro, Brasil.
%   $First version: June 7th, 2003.

if nargin < 5 , U = [1 0 0]; W = [0 1 0]; end
if nargin < 3 , P = [0 0 0]; end
if nargin < 2 , H = 1; end
if nargin < 1 , L = 1; end

% Normalize U:
U = U/norm(U);

% Orthogonalize and then normalize W:
UW = sum(U.*W);
if UW == norm(W);
    disp('Error in RECT: vectors are not linearly independent.');
    return;
end
W = W - UW*U;
W = W/norm(W);

iL = ones(2,1)*L*[-1:2:1]/2;
iH = H*[-1:2:1]'/2*ones(1,2);

x = P(1) + U(1)*iL + W(1)*iH;
y = P(2) + U(2)*iL + W(2)*iH;
z = P(3) + U(3)*iL + W(3)*iH;

if nargout == 0
   surf(x,y,z)
else
   xx = x; yy = y; zz = z;
end
