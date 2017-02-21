function x=GENRANDN(m,R,M)
% this function draws M samples from N(m,R)
% where m is the mean vector(row) and R is the covariance matrix which must
% be positive definite
n=length(m);             % get the dimension
C=chol(R);               % perform cholesky decomp R = C'C
%X=randn(M,n);           % draw M samples from N(0,I)
%x=(X*C+ones(M,1)*m')';
X=randn(n,M);
x=C'*X+m*ones(1,M);