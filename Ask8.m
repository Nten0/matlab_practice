size = [500 1000 5000 10^4 2*10^4];
results = zeros(5,7);

for i = 1:length(size)-3
    n = size(i);
    results(i,1) = n;
    x = ones(n,1);
    A = matrix_build_5853(n,n);
    b = A*x;
    
    x_hat = myCholesky(A,b);
    results(i,2) = norm(b - A*x_hat,2);
    f = @()myCholesky(A,b);    
    results(i,3) = timeit(f);
    
    [x_hat,flag,relres,iter(i,1)] = my_pcg(b,10^(-8),1000,n,n);    
    tmp = matrix_mv_5853(n,x_hat);
    results(i,4) = norm(b - tmp,2);
    f = @()my_pcg(b,10^(-8),1000,n,n);  
    results(i,5) = timeit(f);

    x_hat = calc(A,b);
    results(i,6) = norm(b - A*x_hat,2);
    f = @()calc(A,b);    
    results(i,7) = timeit(f);
end