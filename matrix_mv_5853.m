function [b] = matrix_mv_5853( a,x )
n = length(x);
b = zeros(n,1);

if length(a) == 1
    a = a*ones(n,1);
end

for j = 1:n      
    if (j == 1 || j == n)
        for i = 1:n   
            if (i == 1 && j == 1) || (i == n && j == n)
                b(j) = b(j) + a(i)*x(i);          
            else
                b(j) = b(j) + x(i);
            end
        end
    else
        b(j) = x(1) + x(n) + x(j)*a(j);
    end
end 
end