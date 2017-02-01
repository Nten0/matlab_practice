n = 1000;
A = matrix_build_5853(n,n);
condition = condest(A);
s_eig = sort(eig(A));
histogram(s_eig);