function [xh_k, Ph_k] = ukf_fl(sys, obs, y_k, x_km1, P_km1, Q_k, R_k)
%% UKF
%% Coded by:
% Flavio Eler de Melo (flavio.de-Melo@liverpool.ac.uk)
% University of Liverpool, September 26, 2013
% Usage: [xh_k, Ph_k] = ukf(sys, obs, y_k, x_km1, P_km1, Q_k, R_k)
%
% 'sys' (system state transition) function must take arguments as
% {predicted state vector with added noise} = sys({previous state vector},{multivariate process noise})
% Example: X_k = sys(X_km1, W_k)
%
% 'obs' (observation) functions must take arguments as
% {measurement vector with noise} = obs({state vector},{multivariate measurement noise})
% Example: Y_k = obs(X_k_km1, V_k)
%%
nx = size(x_km1,1);
nu = nx;
ny = size(y_k,1);
nv = ny;

% Parameters for the UKF
na = nx + nu + nv;
alpha = 0.5;
beta = 2;
kappa = 0;
lambda = (alpha^2)*(na + kappa) - na;

% Weights
Wm = zeros(1,2*na+1);
Wc = zeros(1,2*na+1);
Wm(1,1) = lambda/(na + lambda);
Wc(1,1) = Wm(1,1) + (1 -alpha^2 + beta);
Wm(1,2:2*na+1) = 1/(2*(na + lambda));
Wc(1,2:2*na+1) = Wm(1,2:2*na+1);

% Allocate memory
Xa_km1 = zeros(na,2*na+1);
X_k_km1 = zeros(nx,2*na+1);
Y_k_km1 = zeros(ny,2*na+1);
P_k_km1 = zeros(nx,nx);
Pyy = zeros(ny,ny);
Pxy = zeros(nx,ny);

% UKF augmented states
% if det(P_km1) == 0
%     P_km1 = eye(nx)*1e-6;
% end
P_km1
%P_km1 = nearestSPD(P_km1);

Pa_km1 = blkdiag(P_km1, Q_k, R_k)
Xa_km1_m = [x_km1; zeros(nu,1); zeros(nv,1)];



% Scaling parameters and sigma-points
[Si,flag] = chol((na + lambda)*Pa_km1, 'lower');
if flag ~= 0
   SP = chol(P_km1)';
   SQ1 = chol(Q_k(1:3,1:3))';
   %SQ2 = chol(Q_k(4:6,4:6))';
   %SQ3 = chol(Q_k(7:9,7:9))';
   % SQ4 = chol(Q_k(10:12,10:12))';
   % SQ5 = chol(Q_k(13:15,13:15))';
   % SQ6 = chol(Q_k(16:18,16:18))';
    SQ = blkdiag(SQ1);%,SQ2,SQ3,SQ4,SQ5,SQ6
    SR = chol(R_k)';
    Si = blkdiag(SP,SQ,SR);
end
Si
Xa_km1(:,1) = Xa_km1_m;
Xa_km1(:,2:na+1) = Xa_km1_m*ones(1,na) + Si(:,1:na);
Xa_km1(:,na+2:2*na+1) = Xa_km1_m*ones(1,na) - Si(:,1:na);
Xa_km1

% Prediction (time update)
for j = 1:2*na+1
    X_k_km1(1:nx,j) = sys(Xa_km1(1:nx,j), Xa_km1(nx+1:nx+nu,j));
    Y_k_km1(1:ny,j) = obs(X_k_km1(1:nx,j), Xa_km1(nx+nu+1:na,j));
end

x_k_km1 = sum((ones(nx,1)*Wm).*X_k_km1,2);
y_k_km1 = sum((ones(ny,1)*Wm).*Y_k_km1,2);

% Measurement update
for j = 1:2*na+1
    P_k_km1 = P_k_km1 + Wc(1,j)*(X_k_km1(1:nx,j) -x_k_km1)*(X_k_km1(1:nx,j) -x_k_km1)';
    Pyy = Pyy + Wc(1,j)*(Y_k_km1(1:ny,j) -y_k_km1)*(Y_k_km1(1:ny,j) -y_k_km1)';
    Pxy = Pxy + Wc(1,j)*(X_k_km1(1:nx,j) -x_k_km1)*(Y_k_km1(1:ny,j) -y_k_km1)';
end
K_k = Pxy/Pyy;
m_k = x_k_km1 + K_k*(y_k -y_k_km1);
P_k = P_k_km1 - K_k*Pyy*K_k';

% Return filtered state and covariance
xh_k = m_k;
Ph_k = P_k;

return; % bye, bye!!!
