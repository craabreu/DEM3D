function [x,y,z] = particle( N, R, P, u, w, n )

theta = 2*pi*(0:N)/N;
phi   = pi*((0:N)'/N-0.5);
cosphi = cos(phi); cosphi(1) = 0; cosphi(N+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(N+1) = 0;

R_cos_cos = R*cosphi*cos(theta);
R_cos_sin = R*cosphi*sintheta;
R_sin = R*sin(phi)*ones(1,N+1);

x = P(1) + u(1)*R_cos_cos + w(1)*R_cos_sin + n(1)*R_sin;
y = P(2) + u(2)*R_cos_cos + w(2)*R_cos_sin + n(2)*R_sin;
z = P(3) + u(3)*R_cos_cos + w(3)*R_cos_sin + n(3)*R_sin;
