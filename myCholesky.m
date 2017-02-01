function [x] = myCholesky(A,b )
 R = chol(A);
 spy(R)
 x = R\(R'\b);
end