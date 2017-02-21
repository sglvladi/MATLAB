function X=sigmas(x,P,c)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points

A = c*sqrtm(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A]; 