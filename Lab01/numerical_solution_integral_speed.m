clf;

r = 10;
q = 10;
N = 10;
K = 10;

delta_x = 1/N;
delta_y = 1/N;

mu = 0.2;
beta = 0.2;
R = (beta*delta_x)/(2*mu);
%theta = coth(R)-1/R;

theta = 1;
%gamma = R*coth(R);
gamma = 1 + R*theta;

% the grids
x = zeros(1, N+1);
y = zeros(1, K+1);

Psi_x = zeros(N+1, K+1);
Psi_y = zeros(N+1, K+1);
Psi = zeros(N+1, K+1);

function f_val = F(x, y)
  f_val = 0;
endfunction;





for j = 1 : K;
  Psi_x(1,j)=(1+R.*(theta+1)).*(Psi(2,j)-Psi(1,j))./delta_x-delta_x.*(1+theta).*F(x(1,1), y(1,j))/2;
  Psi_x(N,j)=(1+R.*(theta-1)).*(Psi(N-1,j)-Psi(N,j))./delta_x+delta_x.*(1-theta).*F(x(1,N), y(1,j))/2;
endfor;

for i= 2 : N-1;
  for j= 1 : K;
    Psi_x(i,j)=(1-theta).*(1+R.*(theta+1)).*(Psi(i+1,j)-Psi(i,j))./(2*delta_x) + (1+theta).*(1+R.*(theta-1)).*(Psi(i-1,j)-Psi(i,j))./(2*delta_x);
  endfor;
endfor;

for i = 1 : N;
  Psi_y(i,1)=(Psi(i,2)-Psi(i,1))./delta_y - delta_y.*F(x(1,i), y(1,1))./(2*mu);
  Psi_y(i,K)=(Psi(i,K)-Psi(i,K-1))./delta_y + delta_y.*F(x(1,i), y(1,K))./(2*mu);
endfor;

for i = 1 : N;
  for j = 2 : K-1;
    Psi_y(i,j)=(Psi(i,j+1)-Psi(i,j-1))./(2*delta_y);
  endfor;
endfor;

%x = linspace(0, r, 15);
%y = linspace(0, q, 15);
%[xx, yy] = meshgrid(x, y);

%h = quiver(xx, yy, Psi_x, Psi_y, 0.8);

%plot(0,0)
