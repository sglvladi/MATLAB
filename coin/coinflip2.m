numTrials = 20;

x = [1e1 5e1 1e2 5e2 1e3 5e3 1e4 1e5 1e6];
N = numel(x);
heads = zeros(numTrials,N);
tails = heads;

for k = 1:numTrials
    for i=1:N
      for j =1:x(i)
         r = rand;
         if rand <= 0.5
             heads(k,i) = heads(k,i) + 1;
             tails(k,i) = x(i) - heads(k,i);
         end
     end
    end
end

diff = heads-tails;
ratio = heads./tails;

subplot(2,3,1)
semilogx(x,mean(diff),'*-')
title('Mean(Difference)')

subplot(2,3,2)
semilogx(x,std(diff),'*-')
title('StDev(Difference)')

subplot(2,3,3)
semilogx(x,std(diff)./mean(diff),'*-')
title('COV(Difference)')

subplot(2,3,4)
semilogx(x,mean(ratio),'*-')
title('Mean(Ratio)')

subplot(2,3,5)
semilogx(x,std(ratio),'*-')
title('StDev(Ratio)')

subplot(2,3,6)
semilogx(x,std(ratio)./mean(ratio),'*-')
title('COV(Ratio)')

figure

semilogx(x,mean(ratio),'*-')
hold
semilogx(x,mean(ratio)+std(ratio),'r*--')
semilogx(x,mean(ratio)-std(ratio),'r*--')
