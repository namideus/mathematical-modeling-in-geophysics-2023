clf;
F1 = 1;
F2 = 3;
G1 = 3;
G2 = 3;
l = 1;
s = 3;
m = 3;
p = 1;
r = 10;
q = 10;
r_l = (l.*pi)/r;
r_s = (s.*pi)/r;
q_m = (m.*pi)/q;
q_p = (p.*pi)/q;
x = linspace(0, r, 20);
y = linspace(0, q, 20);
[xx, yy] = meshgrid(x, y);
tau_x = (F1.*cos(r_l.*xx) + F2.*sin(r_l.*xx)).*cos(q_m.*yy); % x-component
tau_y = (G1.*cos(r_s.*xx) + G2.*sin(r_s.*xx)).*sin(q_p.*yy); % y-component
h = quiver(xx, yy, tau_x, tau_y, 1.05, 'red');
title ("Wind");
