% gets the OX distance between point p and ELL
%
function [OXdist] = getOX(p,el)

a = el.a;
b = el.b;
C = el.C;
phi = el.phi;
theta = phi*pi/180;
p1= p-C;

rot = [cos(theta) -sin(theta); sin(theta) cos(theta)];
Xrot=rot*p1'; 
x0 = Xrot(1);
y0 = Xrot(2);
c = sqrt((a^2*y0^2) + (b^2*x0^2));

x1 = a*b*x0 / c;
y1 = a*b*y0 / c;

OXdist = norm([x1 y1]);