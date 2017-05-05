function task4(N1,N2,N3,N4)
%Determine the input of the function
N=[N1,N2,N3,N4];
%An array that store the 4 variables in the experiment
P=[0,0,0,0];
%Defualt value of P
%A FOR loop to call function in task 2 with 4 variables
for i=1:4
 P(i)=task2(N(i),5);
 hold on
end
%Change the variables in log form for better display
plot(log10(N),P,'*-');
%Labels's names
xlabel('log10(N)');
ylabel('P');
%Figure title
title('Probabilities P against N in Task4')