%MAE 548 Assignment 6 - Q2
clc, clear
% Step 1 - State the function

%syms d;
%syms v;
%g = 2.65216*(10^(-8))*(d^2.4299)*(v^3.0116) - 0.5;

D = 107.3333;
V = 6.855;
C = [1 0; 0 1]
L = chol(C,"lower")

% Step 2 - Equivalent Transformation

Beta_new = 1.000;
Beta_old = 0.000;
Beta_diff = 1.000;

D_mean  = D
D_std_dev = 10.29023
kesi_D = sqrt(log(1+((D_std_dev/D_mean)^2)))
lamda_D = log(D_mean)-(0.5*(kesi_D^2))

V_mean  = V
V_std_dev = 0.64773
kesi_V = sqrt(log(1+((V_std_dev/V_mean)^2)))
lamda_V = log(V_mean)-(0.5*(kesi_V^2))

iteration = 0;

while abs(Beta_diff) >= 0.0000001

iteration = iteration + 1

D_eqv = D_mean*(1-log(D_mean)+lamda_D)
D_sigma_eqv = D_mean*kesi_D

V_eqv = V_mean*(1-log(V_mean)+lamda_V)
V_sigma_eqv = V_mean*kesi_V

G_mean = 2.65*(10^(-8))*(D_mean^2.4299)*(V_mean^3.0116) - 0.5

% Correlated Parameters

D_corr = (D_mean - D_eqv)/D_sigma_eqv
V_corr = (V_mean - V_eqv)/V_sigma_eqv

% Step 3 - Uncorrelated/Decoupled Reduced Coordinate
%%%% P is a matrix of parameters in the sequence Y,Z,M %%%%

P_prime_prime = L\[D_corr; V_corr]

% Step 4 - Calculate Derivatives

%dg_dd = double(subs(diff(g,d), {d, v}, {D_mean, V_mean}));
%dg_dv = double(subs(diff(g,v), {d, v}, {D_mean, V_mean}));
dg_dd = 2.65216*(10^(-8))*2.4299*(D_mean^1.4299)*(V_mean^3.0116)
dg_dv = 2.65216*(10^(-8))*(D_mean^2.4299)*3.0116*(V_mean^2.0116)

dP_prime_prime = L'*[dg_dd*D_sigma_eqv; dg_dv*V_sigma_eqv]

% Calculation of Gradient Vector

b = 0;
for i = 1:2
    a = (dP_prime_prime(i,1))^2;
    b = b + a;
end
DG2 = b

% Calculation of alpha - direct cosine
DG1 = sqrt(DG2)

for r = 1:2
    alpha(r,1) = dP_prime_prime(r,1)/DG1;
end
alpha

% Step 5 - Calculate New Reduced Coordinate

d = 0;
for j = 1:2
    c = (dP_prime_prime(j,1))*(P_prime_prime(j,1));
    d = (d + c);
end
mul_fac = (d - G_mean)/DG2

for k = 1:2
    P_prime_prime_new(k,1) = mul_fac*dP_prime_prime(k,1);
end
P_prime_prime_new

% Step 6 - Calculation of Beta

f = 0;
for l = 1:2
    e = (P_prime_prime_new(l,1))^2;
    f = f + e;
end

Beta_old = Beta_new
Beta_new = sqrt(f)

% Transfer back to Original Space

P_prime_new = L*P_prime_prime_new

% Correlated to Original

D_new = (P_prime_new(1,1))*D_sigma_eqv + D_eqv
V_new = (P_prime_new(2,1))*V_sigma_eqv + V_eqv

D_mean = D_new
V_mean = V_new

% Step 8 - Checking Beta difference

Beta_diff = Beta_old - Beta_new

% Calculation of Failure Probability
Phi = @(beta) (1- erf(beta/sqrt(2)))/2;
Pf = Phi(Beta_new)

fprintf("========================================================")
end