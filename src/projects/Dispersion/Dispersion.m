clf

% % ASCII
%Vlasiator=load('LDZ2_QPerp_Bxt.dat');
%Real=Vlasiator(:,:);
%Matsize=size(Real);
%lines=Matsize(1);
%cols=Matsize(2);

% % Binary
file = fopen('Dispersion/LDZ2_QPerp_rhot.bin');
Vlasiator = fread(file, 'double');
fclose(file);
% Using the last two numbers for format
lines = Vlasiator(length(Vlasiator) - 1) -1;
cols = Vlasiator(length(Vlasiator));
% Dropping the last two numbers, putting the stream into array shape
Real = reshape(Vlasiator(cols+1:length(Vlasiator)-2), cols, lines);
% Truncating the first line away
Real = Real';

window = hamming(lines)';

for i=1:cols
    Real(:,i) = Real(:, i) .* window;
end


% Partial plotting
k_start = 2;
k_end = 500;
w_start = 1;
w_end = 500;

% Parameters 
length=2.5e8;
ncells=cols;
dx=length/ncells;
dt=0.04;
Temperature=1.0e5;
B0 = 1.001249e-9;
density = 1.0e3;

% Constants
c = 299792458.0;
kb = 1.3806505e-23;
mp = 1.67262171e-27;
me = 9.1093826e-31;
q = 1.60217653e-19;
mu0 = 4*pi()*1.0e-7;
epsilon0 = 8.85418781762e-12;
gamma = 5.0/3.0;

v_th = sqrt(2.0 * kb * Temperature / mp);
w_ci = q*B0/mp;
r_Larmor = v_th / w_ci;


Fourier=fft2(Real);

dk=2*pi() / (cols * dx);
kaxis=0:dk:cols*dk;

dw=2*pi() / (lines * dt);
waxis=0:dw:lines*dw;

imagesc(kaxis(k_start:k_end) * r_Larmor, waxis(w_start:w_end) / w_ci, log10(abs(Fourier(w_start:w_end, k_start:k_end))));
colorbar();
set(gca,'YDir','normal');
hold on

text('Interpreter', 'tex');
xlabel('{\itk}{\cdot}{\itr}_L');
ylabel('\omega/\omega_{ci}');

% Alfvén wave
vA = B0 / sqrt(mu0*density*mp);
%kaxis2 = kaxis.*kaxis;
%kaxis4 = kaxis2.*kaxis2;
%de2 = 5314^2;
%one = ones(1, cols + 1, 'double');
%omega2Val = 0.5 * (kaxis2 ./ (one + kaxis2 * de2) .* (2.0*one + kaxis2 ./ (one + kaxis2 * de2)) - sqrt(kaxis4 ./ ((one + kaxis2*de2).*(one + kaxis2*de2).*(one + kaxis2*de2)) .* (4.0*kaxis2 + kaxis4 ./ (one + kaxis2 * de2))));
plot(kaxis(k_start:k_end) * r_Larmor, vA * kaxis(k_start:k_end));% / sqrt(1.0 + vA*vA / (c*c)))

% Ion-acoustic wave
cS = sqrt(gamma * kb * Temperature / mp);
%Ldebye2 = epsilon0 * kb * Temperature / (density * q * q);
plot(kaxis(k_start:k_end) * r_Larmor, kaxis(k_start:k_end)*cS);% ./ sqrt(1.0 + kaxis.*kaxis*Ldebye2));
%plot(kaxis(k_start:k_end) * r_Larmor, kaxis(k_start:k_end)* sqrt(kb * Temperature / mp)/ w_ci);

% Magnetosonic wave
plot(kaxis(k_start:k_end) * r_Larmor, kaxis(k_start:k_end) * sqrt((cS*cS + vA * vA)));% / (1.0 + vA * vA / (c*c))))
%plot(kaxis, kaxis * 2 * dw * sqrt((cS*cS + vA * vA)));
%plot(kaxis, kaxis * 3 * dw * sqrt((cS*cS + vA * vA)));
%plot(kaxis, kaxis * 4 * dw * sqrt((cS*cS + vA * vA)));

% Numerical propagation
V = dx/dt;
plot(kaxis(k_start:k_end) * r_Larmor, kaxis(k_start:k_end) * V);

% "Langmuir" wave

% Light
%plot(kaxis, kaxis * dw * c)

% Ion cyclotron frequency
plot(kaxis(k_start:k_end) * r_Larmor, 1);
plot(kaxis(k_start:k_end) * r_Larmor, 2);
plot(kaxis(k_start:k_end) * r_Larmor, 3);
plot(kaxis(k_start:k_end) * r_Larmor, 4);

% Ion lower hybrid frequency
w_ce = q*B0/me;
w_pi = sqrt(density * q^2 / (mp * epsilon0));
w_pe = sqrt(density * q^2 / (me * epsilon0));
w_lh = sqrt((w_pi^2 + w_ci^2) / (1 + w_pe^2 / w_ce^2));
%plot(kaxis(k_start:k_end) * r_Larmor, w_lh / w_ci);

% Ion plasma frequency
%plot(kaxis(k_start:k_end), w_pi);

