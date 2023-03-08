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
global r_l = (l.*pi)/r;
global r_s = (s.*pi)/r;
global q_m = (m.*pi)/q;
global q_p = (p.*pi)/q;
global mu = 0.2;
global beta = 0.2;
global N = 50;
global K = 50;
global M = 100;
global delta_x = 1/(N-1);
global delta_y = 1/(K-1);
global R = (beta*delta_x)/(2*mu);
%theta = coth(R)-1/R;
global theta = 1;
global gamma = 1;
%gamma = R*coth(R);
%gamma = 1 + abs(R);

%------------------------------------------------------------------------------------------------------------------------
% Coefficients of the analytical solution
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
  global F1;
  global F2;
  global G1;
  global G2;
  global q_p;
  global r_s;
  global q_m;
  global r_l;
  f_val = -q_m*(F1*cos(r_l*x)+F2*sin(r_l*x))*sin(q_m*y)-r_s*(-G2*cos(r_s*x)+G1*sin(r_s*x)*sin(q_p*y));
endfunction;

%------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------
x1 = linspace(0, r, 50);
y1 = linspace(0, q, 50);
[xx1, yy1] = meshgrid(x1, y1);
Psi_analytic = (C_1.*exp(A.*xx1)+C_2.*exp(B.*xx1)+D_1.*cos(r_l.*xx1)+D_2.*sin(r_l.*xx1)).*sin(q_m.*yy1)
+ (C_1_bar.*exp(A_bar.*xx1)+C_2_bar.*exp(B_bar.*xx1)+D_1_bar.*cos(r_s.*xx1) + D_2_bar.*sin(r_s.*xx1)).*sin(q_p.*yy1);

%------------------------------------------------------------------------------------------------------------------------
Psi_x = zeros(N+1, K+1);
Psi_y = zeros(N+1, K+1);
Psi = zeros(N+1, K+1);
Psi_prev = zeros(N+1, K+1);
[xx, yy] = meshgrid(x, y);

for iter = 1 : M;
  disp("iter"), disp(iter)
  for j = 2 : K-1;
    for i = 2 : N-1;
      Psi(i,j) = ((mu*(gamma - R)*Psi(i - 1,j))/(delta_x*delta_x)+(mu*Psi(i,j-1))/(delta_y*delta_y)+ (mu*(gamma + R)*Psi_prev(i + 1,j))/(delta_x*delta_x)+(mu*Psi_prev(i,j+1))/(delta_y*delta_y))/((2*mu*gamma)/(delta_x*delta_x)+(2*mu)/(delta_y*delta_y))-F(x(1,i),y(1,j));
    endfor;
  endfor;

  clf;
  subplot (1, 2, 1)
  contour(xx1, yy1, Psi_analytic, 25);
  title ("Analytical solution (Flow)");

  subplot (1, 2, 2)
  contour(xx, yy, Psi, 25);
  title ("Numerical solution (Flow)");

  Psi_prev = Psi;
  pause(1);

  % Save previous calculation
  %for j = 1 : K;
   % for i = 1 : N;
    %  Psi_prev(i,j) = Psi(i,j);
    %endfor;
  %endfor;

endfor;
%------------------------------------------------------------------------------------------------------------------------



%for j = 1 : K;
%  Psi_x(1,j)=(1+R.*(theta+1)).*(Psi(2,j)-Psi(1,j))./delta_x-delta_x.*(1+theta).*F(x(1,1), y(1,j))/2;
%  Psi_x(N,j)=(1+R.*(theta-1)).*(Psi(N-1,j)-Psi(N,j))./delta_x+delta_x.*(1-theta).*F(x(1,N), y(1,j))/2;
%endfor;

%for i= 2 : N-1;
%  for j= 1 : K;
%    Psi_x(i,j)=(1-theta).*(1+R.*(theta+1)).*(Psi(i+1,j)-Psi(i,j))./(2*delta_x) + (1+theta).*(1+R.*(theta-1)).*(Psi(i-1,j)-Psi(i,j))./(2*delta_x);
%  endfor;
%endfor;

%for i = 1 : N;
%  Psi_y(i,1)=(Psi(i,2)-Psi(i,1))./delta_y - delta_y.*F(x(1,i), y(1,1))./(2*mu);
%  Psi_y(i,K)=(Psi(i,K)-Psi(i,K-1))./delta_y + delta_y.*F(x(1,i), y(1,K))./(2*mu);
%endfor;

%for i = 1 : N;
%  for j = 2 : K-1;
%    Psi_y(i,j)=(Psi(i,j+1)-Psi(i,j-1))./(2*delta_y);
%  endfor;
%endfor;

%x = linspace(0, r, 15);
%y = linspace(0, q, 15);
%[xx, yy] = meshgrid(x, y);

%h = quiver(xx, yy, Psi_x, Psi_y, 0.8);

%plot(0,0)
