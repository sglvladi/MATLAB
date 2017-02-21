N = 1e2;
C = zeros(N,1);
tic
for i = 1:N
    u = log(sqrt(0.1)*randn + 2);
    sy = log(sqrt(0.2)*randn + 10);
    sx = log(sqrt(0.05)*randn + 5);

   C(i) = 100/(2*pi*u*sqrt(sx*sy)) * exp(-0.5*((40^2/sx) + (100/sy)));
end

P = prctile(C,95)
toc
