n = 1000;
A = matrix_build_5853(n,n);
R = chol(A);
spy(R);
p = amd(A);
A = A(p,p);
R = chol(A);
spy(R);