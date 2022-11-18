%MAE 548 Assignment 6 - Q2
clc, clear
% Step 1 - State the function

syms y;
syms z;
syms m;
g = (y^2)*z - m;

Y = 20;
Z = 25;
M = 8000;
C = [1 0 0.3; 0 1 0; 0.3 0 1]
L = chol(C,"lower")

% Step 2 - Equivalent Transformation

Beta_new = 1.000;
Beta_old = 0.000;
Beta_diff = 1.000;

Y_mean  = Y
Y_std_dev = 3
kesi_Y = sqrt(log(1+((Y_std_dev/Y_mean)^2)))
lamda_Y = log(Y_mean)-(0.5*(kesi_Y^2))

Z_mean = Z
Z_std_dev = 2.5

M_mean  = M
M_std_dev = 1000
kesi_M = sqrt(log(1+((M_std_dev/M_mean)^2)))
lamda_M = log(M_mean)-(0.5*(kesi_M^2))

iteration = 0;

while abs(Beta_diff) >= 0.000000001

iteration = iteration + 1

Y_eqv = Y_mean*(1-log(Y_mean)+lamda_Y)
Y_sigma_eqv = Y_mean*kesi_Y

Z_eqv = Z
Z_sigma_eqv = Z_std_dev

M_eqv = M_mean*(1-log(M_mean)+lamda_M)
M_sigma_eqv = M_mean*kesi_M

G_mean = (Y_mean^2)*Z_mean - M_mean

% Correlated Parameters

Y_corr = (Y_mean - Y_eqv)/Y_sigma_eqv
Z_corr = (Z_mean - Z_eqv)/Z_sigma_eqv
M_corr = (M_mean - M_eqv)/M_sigma_eqv

% Step 3 - Uncorrelated/Decoupled Reduced Coordinate
% P is a matrix of parameters in the sequence Y,Z,M

P_prime_prime = L\[Y_corr; Z_corr; M_corr]

% Step 4 - Calculate Derivatives

dg_dy = double(subs(diff(g,y), {y, z, m}, {Y_mean, Z_mean, M_mean}));
dg_dz = double(subs(diff(g,z), {y, z, m}, {Y_mean, Z_mean, M_mean}));
dg_dm = double(subs(diff(g,m), {y, z, m}, {Y_mean, Z_mean, M_mean}));

dP_prime_prime = L'*[dg_dy*Y_sigma_eqv; dg_dz*Z_sigma_eqv; dg_dm*M_sigma_eqv]

% Calculation of Gradient Vector

b = 0;
for i = 1:3
    a = (dP_prime_prime(i,1))^2;
    b = b + a;
end
DG2 = b

% Calculation of alpha - direct cosine
DG1 = sqrt(DG2)

for r = 1:3
    alpha(r,1) = dP_prime_prime(r,1)/DG1;
end
alpha

% Step 5 - Calculate New Reduced Coordinate

d = 0;
for j = 1:3
    c = (dP_prime_prime(j,1))*(P_prime_prime(j,1));
    d = (d + c);
end
mul_fac = (d - G_mean)/DG2

for k = 1:3
    P_prime_prime_new(k,1) = mul_fac*dP_prime_prime(k,1);
end
P_prime_prime_new

% Step 6 - Calculation of Beta

f = 0;
for l = 1:3
    e = (P_prime_prime_new(l,1))^2;
    f = f + e;
end

Beta_old = Beta_new
Beta_new = sqrt(f)

% Transfer back to Original Space

P_prime_new = L*P_prime_prime_new

% Correlated to Original

Y_new = (P_prime_new(1,1))*Y_sigma_eqv + Y_eqv
Z_new = (P_prime_new(2,1))*Z_sigma_eqv + Z_eqv
M_new = (P_prime_new(3,1))*M_sigma_eqv + M_eqv

Y_mean = Y_new
Z_mean = Z_new
M_mean = M_new

% Step 8 - Checking Beta difference

Beta_diff = Beta_old - Beta_new

% Calculation of Failure Probability
Phi = @(beta) (1- erf(beta/sqrt(2)))/2;
Pf = Phi(Beta_new)

fprintf("========================================================")
end