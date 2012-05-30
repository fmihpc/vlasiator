// Rankine-Hugoniot conditions to determine the right state from the left state of an MHD oblique shock

clf
clear

function result = fct(x)
    result = (V12 - x * vA12)^2 * ...
    (x * vs12 + 0.5 * V12 * cos(theta)^2 * (x * (Gamma - 1) -(Gamma + 1))) + ...
    0.5 * vA12 * V12 * sin(theta)^2 * x * ((Gamma + x * (2 - Gamma)) * V12 - x * vA12 * ((Gamma + 1) - x * (Gamma -1)))
endfunction

// Constants
mp = 1.67e-27;
mu0 = 4.0 * %pi * 1e-7;
kb = 1.3806505e-23;

Vx1 = 1.0e5
Vy1 = 1.0e4
Bx1 = 1.0e-9
By1 = 1.0e-10
T1 = 1.0e5
rho1 = 1e6
Gamma = 5./3.;

theta = acos(Bx1/sqrt(Bx1^2 + By1^2));
vA12 = Bx1^2 / (mp * rho1 * mu0);
V12 = Vx1^2 + Vy1^2;
P1 = rho1 * kb * T1;
vs12 = Gamma * kb * T1 / mp;

X = fsolve(0.0, fct)

rho2 = rho1 * X
Vx2 = Vx1 / X
Vy2 = Vy1 * (V12 - vA12) / (V12 - X * vA12)

V22 = Vx2^2 + Vy2^2;

Bx2 = Bx1
By2 = By1 * (V12 - vA12) * X / (V12 - X * vA12)
P2 = P1 * (X + (Gamma - 1) * X * V12 * (1 - V22 / V12) / (2.0 * vs12));
T2 = P2 / (rho2 * kb)


x=-0:0.1:30;

for i = 1:length(x)
    y(i) = fct(x(i));
end

plot(x, y)
