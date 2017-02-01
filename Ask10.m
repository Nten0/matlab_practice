size = [500 1000 5000 10^4 2*10^4];
results2 = zeros(5,3);

for i = 1:length(size)
    n = size(i);
    results2(i,1) = n;
    A = matrix_build_5853(n,n);
    x = ones(n,1);
    b = matrix_mv_5853(n,x);
    p = amd(A);
    B = A(p,p);
    b_hat = b(p);
    y = myCholesky(B,b_hat);
    x = y(p);
    
    results2(i,2) = norm(b - A*x,2);
    f = @()myCholesky(B,b_hat);    
    results2(i,3) = timeit(f);
end
