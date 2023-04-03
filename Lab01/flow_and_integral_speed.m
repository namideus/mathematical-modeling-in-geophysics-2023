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
mu = 0.2;
beta = 0.2;

A = -beta/(2.*mu)+sqrt((beta/(2.*mu)).^2+q_m.^2);
B = -beta/(2.*mu)-sqrt((beta/(2.*mu)).^2+q_m.^2);
D_1 = q_m.*(mu.*(r_l.^2+q_m.^2)*F1 + beta.*r_l.*F2)/(mu.^2.*(r_l.^2+q_m.^2).^2+beta.^2.*r_l.^2);
D_2 = q_m.*(mu.*(r_l.^2+q_m.^2)*F2 - beta.*r_l.*F1)/(mu.^2.*(r_l.^2+q_m.^2).^2+beta.^2.*r_l.^2);
C_1 = D_1*(exp(B.*r)-(-1).^l)/(exp(A.*r)-exp(B.*r));
C_2 = D_1*((-1).^l-exp(A.*r))/(exp(A.*r)-exp(B.*r));

A_bar = -beta/(2.*mu)+sqrt((beta/(2.*mu)).^2+q_p.^2);
B_bar = -beta/(2.*mu)-sqrt((beta/(2.*mu)).^2+q_p.^2);
D_1_bar = r_s.*(mu.*(r_s.^2+q_p.^2)*G2 - beta.*r_s.*G1)/(mu.^2.*(r_s.^2+q_p.^2).^2+beta.^2.*r_s.^2);
D_2_bar = -r_s.*(mu.*(r_s.^2+q_p.^2)*G1 + beta.*r_s.*G2)/(mu.^2.*(r_s.^2+q_p.^2).^2+beta.^2.*r_s.^2);
C_1_bar = D_1_bar*(exp(B_bar.*r)-(-1).^s)/(exp(A_bar.*r)-exp(B_bar.*r));
C_2_bar = D_1_bar*((-1).^s-exp(A_bar.*r))/(exp(A_bar.*r)-exp(B_bar.*r));

x = linspace(0, r, 50);
y = linspace(0, q, 50);
[xx, yy] = meshgrid(x, y);

Psi = (C_1.*exp(A.*xx)+C_2.*exp(B.*xx)+D_1.*cos(r_l.*xx)+D_2.*sin(r_l.*xx)).*sin(q_m.*yy)+ (C_1_bar.*exp(A_bar.*xx)+C_2_bar.*exp(B_bar.*xx)+D_1_bar.*cos(r_s.*xx) + D_2_bar.*sin(r_s.*xx)).*sin(q_p.*yy);

subplot (1, 2, 1)
contour(xx, yy, Psi, 25);
title ("Flow");

x1 = linspace(0, r, 20);
y1 = linspace(0, q, 20);
[xx, yy] = meshgrid(x1, y1);

%[xx, yy] = meshgrid(0:0.5:r, 0:0.5:q);
% v-component
V=(A.*C_1.*exp(A.*xx)+B.*C_2.*exp(B.*xx)-r_l.*D_1.*sin(r_l.*xx)+r_l.*D_2.*cos(r_l.*xx)).*sin(q_m.*yy)+(A_bar.*C_1_bar.*exp(A_bar.*xx)+B_bar.*C_2_bar.*exp(B_bar.*xx)-r_s.*D_1_bar.*sin(r_s.*xx)+r_s.*D_2_bar.*cos(r_s.*xx)).*sin(q_p.*yy);
V=V.*(-1);
% u-component
U=q_m.*(C_1.*exp(A.*xx)+C_2.*exp(B.*xx)+D_1.*cos(r_l.*xx)+D_2.*sin(r_l.*xx)).*cos(q_m.*yy)+q_p.*(C_1_bar.*exp(A_bar.*xx)+C_2_bar.*exp(B_bar.*xx)+D_1_bar.*cos(r_s.*xx) + D_2_bar.*sin(r_s.*xx)).*cos(q_p.*yy);
subplot (1, 2, 2)
h = quiver(xx, yy, U, V, 1);
title ("Integral speed");
