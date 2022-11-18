%MAE 548 Assignment 6 - Q3
%clc, clear

Y_mean = 20;
Y_sigma = 3;
Z_mean = 25;
Z_sigma = 2.5;
M_mean = 8000;
M_sigma = 1000;
kesi_Y = sqrt(log(((Y_sigma/Y_mean)^2)+1));
kesi_M = sqrt(log(((M_sigma/M_mean)^2)+1));
lamda_Y = log(Y_mean) - 0.5*(kesi_Y^2);
lamda_M = log(M_mean) - 0.5*(kesi_M^2);

n = 100000;
k = zeros(1,n);

for i = 1: n
    Y = lognrnd (lamda_Y, kesi_Y);
    Z = normrnd (25, 2.5);
    M = lognrnd (lamda_M, kesi_M);
    G = Y^2*Z - M;
    if -1000 <= Y^2*Z - M & 1000 >= Y^2*Z - M
        k(i) = 1;
    end
end

s = sum(k);
p = s/n

fprintf('===============================================')