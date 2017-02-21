function [m] = throwNeedles2(numNeedles)

inCircle = 0;
P = lhsdesign(numNeedles,2);
for i = 1:numNeedles
    if sqrt(P(i,1)^2+P(i,2)^2) <= 1;
        inCircle = inCircle + 1;
    end
end
m = 4*(inCircle/numNeedles);
