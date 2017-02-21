x = [1e1 5e1 1e2 5e2 1e3 5e3 1e4 1e5 1e6];
N = numel(x);
heads = zeros(N,1);
tails = heads;

for i=1:N
  for j =1:x(i)
     r = rand;
     if rand <= 0.5
         heads(i) = heads(i) + 1;
         tails(i) = x(i) - heads(i);
     end
 end
end

subplot(2,1,1)
semilogx(x,heads-tails,'*-')
title('Difference')
subplot(2,1,2)
semilogx(x,heads./tails,'*-')
title('Ratio')


