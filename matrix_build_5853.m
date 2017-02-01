function [A] = matrix_build_5853(n,a)
A = sparse(n,n);
for i = 2:n
    A(1,i) = 1;
    A(n,i-1) = 1;
    A(i,1) = 1;
    A(i-1,n) = 1;
end

if length(a) == 1
    for i = 1:n
        A(i,i) = a;
    end
else
     for i = 1:n
        A(i,i) = a(i);
    end
end
end