clf;

%------------------------------------------------------------------------------------------------------------------------
global F1 = 1;
global F2 = 3;
global G1 = 3;
global G2 = 3;
global l = 1;
global s = 3;
global m = 3;
global p = 1;
global r = 10;
global q = 10;
global r_l = (l*pi)/r;
global r_s = (s*pi)/r;
global q_m = (m*pi)/q;
global q_p = (p*pi)/q;
global mu = 0.2;
global beta = 0.2;
global N = 64;
global K = 64;
global M = 1500;
global delta_x = r/(N-1);
global delta_y = q/(K-1);
global R = (beta*delta_x)/(2*mu);
%global theta = 1;
global theta = coth(R)-1/R;
%global gamma = 1;
global gamma = R*coth(R);
%gamma = 1 + abs(R);

%------------------------------------------------------------------------------------------------------------------------
% Coefficients of the analytical solution
global A = -beta/(2.*mu)+sqrt((beta/(2.*mu)).^2+q_m.^2);
global B = -beta/(2.*mu)-sqrt((beta/(2.*mu)).^2+q_m.^2);
global D_1 = q_m.*(mu.*(r_l.^2+q_m.^2)*F1 + beta.*r_l.*F2)/(mu.^2.*(r_l.^2+q_m.^2).^2+beta.^2.*r_l.^2);
global D_2 = q_m.*(mu.*(r_l.^2+q_m.^2)*F2 - beta.*r_l.*F1)/(mu.^2.*(r_l.^2+q_m.^2).^2+beta.^2.*r_l.^2);
global C_1 = D_1*(exp(B.*r)-(-1).^l)/(exp(A.*r)-exp(B.*r));
global C_2 = D_1*((-1).^l-exp(A.*r))/(exp(A.*r)-exp(B.*r));

global A_bar = -beta/(2.*mu)+sqrt((beta/(2.*mu)).^2+q_p.^2);
global B_bar = -beta/(2.*mu)-sqrt((beta/(2.*mu)).^2+q_p.^2);
global D_1_bar = r_s.*(mu.*(r_s.^2+q_p.^2)*G2 - beta.*r_s.*G1)/(mu.^2.*(r_s.^2+q_p.^2).^2+beta.^2.*r_s.^2);
global D_2_bar = -r_s.*(mu.*(r_s.^2+q_p.^2)*G1 + beta.*r_s.*G2)/(mu.^2.*(r_s.^2+q_p.^2).^2+beta.^2.*r_s.^2);
global C_1_bar = D_1_bar*(exp(B_bar.*r)-(-1).^s)/(exp(A_bar.*r)-exp(B_bar.*r));
global C_2_bar = D_1_bar*((-1).^s-exp(A_bar.*r))/(exp(A_bar.*r)-exp(B_bar.*r));
%------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------
% the grids

x = zeros(1, N+1);
y = zeros(1, K+1);

for i = 1 : N;
  x(1,i)=(i-1)*delta_x;
endfor;

for j = 1 : K;
  y(1,j)=(j-1)*delta_y;
endfor;

function f_val = F(x, y)
  global F1 F2 G1 G2 q_p r_s q_m r_l;
  f_val = -q_m*(F1*cos(r_l*x)+F2*sin(r_l*x))*sin(q_m*y)-r_s*(-G1*sin(r_s*x)+G2*cos(r_s*x))*sin(q_p*y);
endfunction;

function ret = Psi_func(x, y)
  global C_1 A C_2 B D_1 r_l D_2 q_m C_1_bar A_bar C_2_bar B_bar D_1_bar r_s D_2_bar q_p;
  ret=(C_1*exp(A*x)+C_2*exp(B*x)+D_1*cos(r_l*x)+D_2*sin(r_l*x))*sin(q_m*y)+(C_1_bar*exp(A_bar*x)+C_2_bar*exp(B_bar*x)+D_1_bar*cos(r_s*x)+D_2_bar*sin(r_s*x))*sin(q_p*y);
endfunction;

%------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------
x1 = linspace(0, r, 50);
y1 = linspace(0, q, 50);
[xx1, yy1] = meshgrid(x1, y1);
Psi_analytic = (C_1.*exp(A.*xx1)+C_2.*exp(B.*xx1)+D_1.*cos(r_l.*xx1)+D_2.*sin(r_l.*xx1)).*sin(q_m.*yy1)+ (C_1_bar.*exp(A_bar.*xx1)+C_2_bar.*exp(B_bar.*xx1)+D_1_bar.*cos(r_s.*xx1) + D_2_bar.*sin(r_s.*xx1)).*sin(q_p.*yy1);

%------------------------------------------------------------------------------------------------------------------------
Psi = zeros(N+1, K+1);
Psi_prev = zeros(N+1, K+1);
[xx, yy] = meshgrid(x, y);
%------------------------------------------------------------------------------------------------------------------------
for iter = 1 : M;
  disp("iter"), disp(iter)
  for j = 2 : K-1;
    for i = 2 : N-1;
      Psi(i,j)=power(((2*mu*gamma)/(delta_x^2)+(2*mu)/(delta_y^2)),-1)*((mu*(gamma-R)*Psi(i-1,j))/(delta_x^2)+(mu*Psi(i,j-1))/(delta_y^2)+(mu*(gamma+R)*Psi_prev(i+1,j))/(delta_x^2)+(mu*Psi_prev(i,j+1))/(delta_y^2)-F(x(1,i),y(1,j)));
    endfor;
  endfor;

  Psi_prev(1:end,1:end) = Psi(1:end,1:end);

  maxError = 0;
  for j=2:K;
    for i=2:N;
      maxError=max(abs(Psi_func(x(1,i),y(1,j))-Psi(i,j)),maxError);
    endfor;
  endfor;
  disp('Max error: '),disp(maxError);

  subplot (1, 2, 1)
  contour(xx1, yy1, Psi_analytic, 25);
  title ("Analytical solution (Flow)");

  subplot (1, 2, 2)
  contour(yy, xx, Psi, 25);
  title ("Numerical solution (Flow)");
  text (-3.0, -0.7, strcat("max error:\t", num2str(maxError)),"fontsize", 20);

  pause(0.001);
endfor;
%------------------------------------------------------------------------------------------------------------------------

