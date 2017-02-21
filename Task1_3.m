function p = Task1_3(N,R)
% Date Created: 21-10-2016
% Author: Bingkun Guo & Zeen Guo
% Constants:
L=4;
W=2;
r=sqrt((L/2)^2+(W/2)^2);
sum=0;
% Define variables to support plotting
xset_c=zeros(1,N);
yset_c=zeros(1,N);
xset_s=zeros(1,N);
yset_s=zeros(1,N);
for k=1:R
 t=0;
 j=0;
 for i=1:N

 % Use polar coodinates to set random positions
 angle=unifrnd(0,2*pi);
 d=sqrt(unifrnd(0,1))*r;

 % Transfer polar coordinates into rectangular coordinates
 x=d*cos(angle);
 y=d*sin(angle);

 % Check if the positon was goal and put the value into array
 if(abs(x)<(L/2)&&abs(y)<(W/2))
 t=t+1;
 xset_s(1,t)=x;
 yset_s(1,t)=y;
 else
 j=j+1;
 xset_c(1,j)=x;
 yset_c(1,j)=y;
 end;
 end;
 sum = sum +t/N;

end;
% Produce the scatter plot
if(R==1)
 plot(xset_c,yset_c,'bo',xset_s,yset_s,'rx');
 legend('Miss','Score');
 title('Scatter plot(Guo)');
end;
% Calculate the probability
p=sum/R;