clf;

%------------------------------------------------------------------------------------------------------------------------
global f1 = 1;
global f2 = 10;
global omega = 5;
global mu0 = 4*pi*10^(-7);
global sigma = 0.01;
global zmax = 500;
global k =(1-i)*sqrt(omega*mu0*sigma/2);
global N = 500;
global delta_z = zmax/(N-1);
global theta = 1;

maxHrealError = 0;
maxHimageError = 0;
maxErealError = 0;
maxEimageError = 0;
%------------------------------------------------------------------------------------------------------------------------
% the grids
H_num = zeros(N);
E_hum = zeros(N);
z_grid = zeros(N);

for j = 1 : N;
  z_grid(j)=(j-1)*delta_z;
endfor;
%------------------------------------------------------------------------------------------------------------------------
% Exact H(z)
function ret = H(z)
  global k f1 f2 zmax;
  ret=f1*cosh(k*z)-(sinh(k*z)*(f1*cosh(k*zmax)-f2))/sinh(k*zmax);
endfunction;

% Exact E(z)
function ret = E(z)
  global k f1 f2 zmax omega mu0 sigma;
  ret=(f1*sinh(k*z)*(i*omega*mu0))/k+(cosh(k*z)*(f1*cosh(k*zmax)-f2))/((sigma/k)*sinh(k*zmax));
endfunction;
%------------------------------------------------------------------------------------------------------------------------
z = linspace(0, zmax, 50);
%------------------------------------------------------------------------------------------------------------------------
H_num(1) = f1;
E_num(1) = (f1*cosh(k*zmax)-f2)/((sigma/k)*sinh(k*zmax));
H_num(N) = f2;

for j = 1 : N-1;
  H_num(j+1) = H_num(j)-delta_z*sigma*E_num(j);
  E_num(j+1) = E_num(j)+delta_z*i*omega*mu0*H_num(j);
endfor;
%--------------------------------------------extract real and imaginary parts---------------------------------------------
% exact
H_rr = real(H(z));
H_ii = imag(H(z));
E_rr = real(E(z));
E_ii = imag(E(z));

% numerical
Hnum_rr = real(H_num);
Hnum_ii = imag(H_num);
Enum_rr = real(E_num);
Enum_ii = imag(E_num);

% error
maxHreal = max(H_rr);
maxHimage = max(H_ii);
maxEreal = max(E_rr);
maxEimage = max(E_ii);

for j=1:N;
  maxHrealError=max(abs(real(H(z_grid(j)))-Hnum_rr(j))*100/abs(maxHreal),maxHrealError);
  maxHimageError=max(abs(imag(H(z_grid(j)))-Hnum_ii(j))*100/abs(maxHimage),maxHimageError);
  maxErealError=max(abs(real(E(z_grid(j)))-Enum_rr(j))*100/abs(maxEreal),maxErealError);
  maxEimageError=max(abs(imag(E(z_grid(j)))-Enum_ii(j))*100/abs(maxEimage),maxEimageError);
endfor;

% plotting
%  real part of H(z)
subplot (2, 2, 1)
plot(z, H_rr, 'b');
xlabel("z axis");
ylabel("real(H) axis");

hold on;
plot(z_grid, Hnum_rr, 'r');
hold off;
title ("Real part of H(z)");
text (0, -2, strcat("error real(H):\t", num2str(maxHrealError), "%\n", "error imag(H):\t", num2str(maxHimageError), "%\n", ...
"error real(E):\t", num2str(maxErealError), "%\n", "error imag(E):\t", num2str(maxEimageError),"%"),"fontsize", 12);

%  imaginary part of H(z)
subplot (2, 2, 2)
plot(z, H_ii, 'b');
xlabel("z axis");
ylabel("imag(H) axis");

hold on;
plot(z_grid, Hnum_ii, 'r');
hold off;
title ("Imaginary part of H(z)");
text (0, -2, strcat("Error:\t", num2str(100), "%"),"fontsize", 15);

%  real part of E(z)
subplot (2, 2, 3)
plot(z, E_rr, 'b');
xlabel("z axis");
ylabel("real(E) axis");

hold on;
plot(z_grid, Enum_rr, 'r');
hold off;
title ("Real part of E(z)");

%  imaginary part of E(z)
subplot (2, 2, 4)
plot(z, E_ii, 'b');
xlabel("z axis");
ylabel("imag(E) axis");

hold on;
plot(z_grid, Enum_ii, 'r');
hold off;
title ("Imaginary part of E(z)");

disp(maxHrealError);
disp(maxHimageError);
disp(maxErealError);
disp(maxEimageError);
%------------------------------------------------------------------------------------------------------------------------


