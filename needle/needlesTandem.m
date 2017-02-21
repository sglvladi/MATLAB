clear all
close all

N = 1e2;
n = 1e3;
h1 = zeros(N,1);
h2 = h1;

tic
for i = 1:N
    h1(i) = throwNeedles(n,0);
    h2(i) = throwNeedles2(n);
end

subplot(2,1,1)
hist(h1);
title(['Uniform Sampling, Var = ',num2str(var(h1))])
lims = get(gca,'xlim')
subplot(2,1,2)
hist(h2);
title(['LH Sampling, Var = ',num2str(var(h2))])
set(gca,'xlim',lims)
toc