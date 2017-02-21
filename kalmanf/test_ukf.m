clear
n=3;      %number of state
q=0.1;    %std of process 
r=1;    %std of measurement
s.Q=q^2*eye(n); % covariance of process
s.R=r^2*eye(n);        % covariance of measurement  
s.sys=@(x)[sin(x(2));x(2)+1;0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
s.obs=@(x)[x(1);x(2);x(3)];                               % measurement equation
st=[0;0;1];                                % initial state
s.x=st+q*randn(3,1); %initial state          % initial state with noise
s.P = eye(n);                               % initial state covraiance
N=100;                                     % total dynamic steps
xV = zeros(n,N+1);          %estmate        % allocate memory
PV = zeros(1,N+1);
sV = zeros(n,N+1);          %actual
zV = zeros(n,N+1);
ukf = UKalmanFilter(s, 0.5, 0, 2);
for k=1:N
  ukf.s.z = ukf.s.sys(st) + r*randn                     % measurments
  sV(:,k)= st;                             % save actual state
  zV(:,k)  = ukf.s.z;                             % save measurment
  ukf.s = ukf.Iterate(ukf.s);            % ekf 
  xV(:,k+1) = ukf.s.x;                            % save estimate
  PV(k) = ukf.s.P(1,1);
  st = ukf.s.sys(st);                % update process 
end
figure
for k=1:3                                 % plot results
  subplot(4,1,k)
  plot(1:N+1, sV(k,:), '-', 1:N+1, xV(k,:), '--', 1:N+1, zV(k,:), 'r.')
end
subplot(4,1,4)
plot(1:N+1, PV(:,:))
