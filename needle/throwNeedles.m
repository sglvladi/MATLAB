function [m] = throwNeedles(numNeedles,flag)

X = zeros(numNeedles,1);
Y = X;
Xno = X;
Yno = X;

inCircle = 0;
for i = 1:numNeedles
    x = rand;
    y = rand;
    if sqrt(x^2+y^2) <= 1;
        X(i) = x;
        Y(i) = y;
        inCircle = inCircle + 1;
    else
        Xno(i) = x;
        Yno(i) = y;
    end
end
m = 4*(inCircle/numNeedles);

if flag == 1
    hold on
    cx=0:0.025:1;
    cy = sqrt(1-(cx.^2));
    figure(1);
    plot(cx, cy, 'k','linewidth',2);
    title(['Pi estimate = ', num2str(m)])

    plot(X,Y,'.b')
    plot(Xno,Yno,'.r')

    axis square
end


