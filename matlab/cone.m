function [x,y,z] = cone(N, L, Ra, Rp, P, u, w, n, Half)

if nargin < 9
    Half = 0;
end

if Half
    N = N / 2;
    theta  = -pi*(0:N)'/N;
else
    theta  = 2*pi*(0:N)'/N;
end
lambda = L/2*[-1 1];
R = [Ra Rp];

costheta = cos(theta); costheta(1)=1;
sintheta = sin(theta); sintheta(1)=0;
if ~Half
    costheta(N+1)=1;
    sintheta(N+1)=0;
end
lambda   = ones(N+1,1)*lambda;

x = P(1) + u(1)*lambda + (w(1)*costheta + n(1)*sintheta)*R;
y = P(2) + u(2)*lambda + (w(2)*costheta + n(2)*sintheta)*R;
z = P(3) + u(3)*lambda + (w(3)*costheta + n(3)*sintheta)*R;

