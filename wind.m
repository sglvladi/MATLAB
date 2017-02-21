clear all

tic
N = 1e5;

Y = zeros(100,1);
X = zeros(N,1);
E = zeros(N,1);

Y(1) = 1.5;
Y(2) = Y(1);

for i = 1:N
    for t = 3:100
        e = sqrt(0.02)*randn;
        Y(t) = 0.7*Y(t-1) + 0.45*Y(t-2) + e;
    end
    X(i) = sum((Y<log(5.556)));
    if X(i) > 10
        E(i) = 1;
    else
        E(i) = 0;
    end
end
p = sum(E)/N
toc
