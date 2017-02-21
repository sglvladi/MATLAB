N = 1e5;
n = 90;
Z = zeros(N,1);
tic
for i = 1:N
    X = randn(n,1);
    Z(i) = max(X) - min(X);
end
V = var(Z)
toc
