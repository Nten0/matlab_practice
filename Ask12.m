a = 999;
n = 1000;
x = ones(n,1);
b = matrix_mv_5853(n,x);
A = matrix_build_5853(1000,a);
condition = condest(A);
[x_hat,flag,relres,iterations] = my_pcg(b,10^(-8),1000,n,a)