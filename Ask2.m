n = 6;
a = 5853*ones(n,1);
A = matrix_build_5853(n,a);

c = 1;
index = 2;
IA(1) = 1;

for i = 1:length(A)
    for j = 1:length(A)
        if A(i,j) ~= 0
            VAL(c,1) = A(i,j);
            COL(c,1) = j;
            c = c + 1;
        end
    end
    IA(index,1) = c;
    index = index + 1;
end