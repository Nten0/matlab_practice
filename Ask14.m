n = 1000;
A = sparse(n,n);
F = zeros(n,1);

for i = 1:n
    if i == 1
        A(1,1) = 2*(n+1)^2 + 6;
        A(1,2) = -(n+1)^2 - 2;
        F(1) = 1;
    elseif i == n
        A(n,n-1) = -(n+1)^2 - n^2 + n ;
        A(n,n) = 2*(n+1)^2 + 2*n^2 + 4;
        F(n) = 1 + (n+1)^2 + n^2 + n;
    else
        A(i,i-1) = - ((n+1)^2 + i^2 - i);
        A(i,i) = 2*(n+1)^2 + 2*i^2 + 4;
        A(i,i+1) = -(n+1)^2 - i^2 - i;
        F(i) = 1;
    end
end
[X,flag,relres,iterations] = pcg(A,F,[],1000);