n = 6;
a = 5853*ones(n,1);
A = matrix_build_5853(n,a);
B = tril(A);

c = 1;
index = 2;
IA(1) = 1;

for i = 1:length(A)
    for j = 1:length(A)
        if A(j,i) ~= 0
            VAL(c,1) = A(j,i);
            ROW(c,1) = j;
            c = c + 1;
        end
    end
    IA(index,1) = c;
    index = index + 1;
end