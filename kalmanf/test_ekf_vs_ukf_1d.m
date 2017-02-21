%% Clear everything
%clear

%pendulum_sim.m

%% Initiate KF parameters
n=2;      %number of state
Dt=0.01;
q=sqrt(0.02);    %std of process 
r=sqrt(0.1);    %std of measurement
s.Q=[Dt^3/3 Dt^2/2; Dt^2/2 Dt]*q^2; % covariance of process
s.R=r^2*eye(1);        % covariance of measurement  
s.sys=@(t)(@(x)[x(1)+x(2)*t;x(2)-9.81*sin(x(1))*t]);  % nonlinear state equations
s.obs=@(x) sin(x(1));                               % measurement equation
st=[1.6;0];                                % initial state
s.x=st;%+mvnrnd([0;0],s.Q)'; %initial state          % initial state with noise
s.P = eye(2)*r^2;                               % initial state covraiance
x_ukf = s.x;
P_ukf = s.P;
N=900;                                     % total dynamic steps

%% Instantiate UKF and vars to store output
xV_ukf = zeros(n,N);          %estmate        % allocate memory
xV_ukf(:,1) = s.x;
PV_ukf = zeros(1,N);
%PV_ukf = cell(1,N);    % use to display ellipses 
sV_ukf = zeros(n,N);          %actual
sV_ukf(:,1) = [1.6;0]
zV_ukf = zeros(n,N);
eV_ukf = zeros(n,N);
ukf = UKalmanFilter(s, 0.5, 1, 2);

%% Instantiate EKF and vars to store output
xV_ekf = zeros(n,N);          %estmate        % allocate memory
xV_ekf(:,1) = s.x;
PV_ekf = zeros(1,N);    
%PV_ekf = cell(1,N);     % use to display ellipses 
sV_ekf = zeros(n,N);          %actual
zV_ekf = zeros(n,N);
eV_ekf = zeros(n,N);
ekf = EKalmanFilter(s);

%% Initiate PF parameters

% Process equation x[k] = sys(k, x[k-1], u[k]);
nx = 2;  % number of states
sys = @(t)(@(k, xkm1, uk)[xkm1(1)+xkm1(2)*t+uk(1);xkm1(2)-9.81*sin(xkm1(1))*t+uk(2)]); % (returns column vector)

% Observation equation y[k] = obs(k, x[k], v[k]);
ny = 1;                                           % number of observations
obs = @(k, xk, vk) [sin(xk(1))+vk(1)];                  % (returns column vector)

% PDF of process noise and noise generator function
nu = 2;                                           % size of the vector of process noise
sigma_u = sqrt(0.02);
cov_u = [Dt^3/3 Dt^2/2; Dt^2/2 Dt]*sigma_u^2;
p_sys_noise   = @(u) mvnpdf(u, zeros(1, nu), cov_u);
gen_sys_noise = @(u) mvnrnd(zeros(1, nu), cov_u);         % sample from p_sys_noise (returns column vector)

% PDF of observation noise and noise generator function
nv = 1;                                           % size of the vector of observation noise
sigma_v = sqrt(0.1);
cov_v = sigma_v^2*eye(nv);
p_obs_noise   = @(v) mvnpdf(v, zeros(1, nv), cov_v);
gen_obs_noise = @(v) mvnrnd(zeros(1, nv), cov_v);         % sample from p_obs_noise (returns column vector)

% Initial PDF
% p_x0 = @(x) normpdf(x, 0,sqrt(10));             % initial pdf
gen_x0 = @(x) mvnrnd([1.6;0],cov_u);               % sample from p_x0 (returns column vector)

% Transition prior PDF p(x[k] | x[k-1])
% (under the suposition of additive process noise)
% p_xk_given_xkm1 = @(k, xk, xkm1) p_sys_noise(xk - sys(k, xkm1, 0));

% 
% (under the suposition of additive process noise)
p_yk_given_xk = @(k, yk, xk) p_obs_noise((yk - obs(k, xk, zeros(1, nv)))');

% Number of time steps
T = 900;

% Separate memory space
%x = [x_true(:,1), y_true(:,1)]';
%y = [obs_x'; obs_y']; % True state and observations
xh = zeros(nx, T); yh = zeros(ny, T); % Filtered state and observations

% Initial state and observatiopns
%xh0 = [obs_x(2); obs_y(2); obs_x(2)-obs_x(1); obs_y(2)-obs_y(1)] ; % 2-point state initialisation
% Single-point initiation
xh0=[1.6;0]; %initial state
xh(:,1) = xh0; yh(:,1) = obs(1, xh0, zeros(1, nv)); % Filtered

% Assign PF parameter values
pf.k               = 1;                   % initial iteration number
pf.Np              = 1000;                 % number of particles
%pf.w               = zeros(pf.Np, T);     % weights
pf.particles       = zeros(nx, pf.Np, T); % particles
pf.gen_x0          = gen_x0;              % function for sampling from initial pdf p_x0 (used to generate initial particles)
pf.obs             = p_yk_given_xk;       % function of the observation likelihood PDF p(y[k] | x[k])
pf.sys_noise       = gen_sys_noise;       % function for generating system noise
%pf.p_x0 = p_x0;                          % initial prior PDF p(x[0])
%pf.p_xk_given_ xkm1 = p_xk_given_xkm1;   % transition prior PDF p(x[k] | x[k-1])
pf.xhk = xh0;
pf.sys = sys(1);
pf.resampling_strategy = 'systematic_resampling';
my_pf = ParticleFilter(pf);


x = zeros(nx, T);
x(:,1)=st;
eV_obs =[];
for k=2:N
    fprintf('Iteration = %d/%d\n',k,N);
    ukf.s.sys = s.sys(Dt);
    ekf.s.sys = s.sys(Dt);
    my_pf.pf.k = k;
    my_pf.pf.sys = sys(Dt);
    %% Generate new state
    st = X(:,k);
    %% Get next measurement
    ukf.s.z = Y(:,k);                     % measurments
    ekf.s.z = ukf.s.z;
    my_pf.pf.z = ukf.s.z(1);

    %% Store new state and measurement
    sV_ukf(:,k)= st;                             % save actual state
    zV_ukf(:,k)  = ukf.s.z;                             % save measurment
    sV_ekf(:,k)= st;                             % save actual state
    zV_ekf(:,k)  = ekf.s.z; 

    %% Compute Q
    %ukf.s.Q = eye(2)*q^2;
    %ekf.s.Q = ukf.s.Q;
    
    %% Iterate both filters
    ekf.s = ekf.Iterate(ekf.s);
    %ukf.s.sys = s.sys(N);
    %[x_ukf, P_ukf] = ukf_fl(ukf_sys, ukf_obs, ekf.s.z, x_ukf, P_ukf,  s.Q, s.R);
    ukf.s = ukf.Iterate(ukf.s);            % ekf 
    my_pf.pf = my_pf.Iterate(my_pf.pf);
    
    %% Store estimated state and covariance
    xV_ukf(:,k) = ukf.s.x;%x_ukf;                            % save estimate
    PV_ukf(k)= ukf.s.P(1,1); % P_ukf(1,1) ; 
    %PV_ukf{k}= ukf.s.P;    % Use to store whole covariance matrix
    xV_ekf(:,k) = ekf.s.x;                            % save estimate
    PV_ekf(k) = ekf.s.P(1,1);
    %PV_ekf{k} = ekf.s.P;    % Use to store whole covariance matrix
    xh(:,k) = my_pf.pf.xhk(:,k);
    x(:,k)=st;
    %% Compute squared error
    eV_ukf(:,k) = (ukf.s.x - st).*(ukf.s.x - st);
    eV_ekf(:,k) = (ekf.s.x - st).*(ekf.s.x - st);
    eV_obs(:,k) = (ukf.s.z - st).*(ukf.s.z - st);
                    % update process 
    
end

%% Compute & Print RMSE
RMSE_obs = sqrt(sum(eV_obs,2)/N)
RMSE_ukf = sqrt(sum(eV_ukf,2)/N)
RMSE_ekf = sqrt(sum(eV_ekf,2)/N)
%% Compute RMSE
err = (xh(:,:) - x).*(xh(:,:) - x);
RMSE_pf = sqrt(sum(err,2)/T)

%% Plot results
figure
%for k=1:2                                 % plot results
%     subplot(2,3,1)
% %     figure
%      hold on
%     plot(1:N, sV_ukf(1,:), 'g-', 'LineWidth', 3)
%     plot(1:N, zV_ukf(1,:), 'r.', 'LineWidth', 1)
%     plot(1:N, xV_ukf(1,:), 'b-', 'LineWidth', 1)
%     str = sprintf('UKF estimated state vs. Time');
%     title(str)
%     h_legend = legend('Real', 'Meas', 'UKF');
%     set(h_legend,'FontSize',5);
%     axis([0,500,-4,4])

%end
% subplot(2,1,2)
% plot(1:N, eV_ukf(1,:))
% title(sprintf('UKF RMSE vs. Time'))
% legend('UKF');
% axis([0,500,0,0.3])

%% Plot results
%figure
%for k=1:2   
% colourMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
% colourMap('y') = [ 0.9290    0.6940    0.1250]
% plot results
    subplot(2,3,[1:3])
%     figure
     hold on
    plot(Dt*([1:N]), zV_ukf(1,:), 'color', [0.7 0.7 0.7], 'Marker', '.')
    plot(Dt*(1:N), sV_ukf(1,:), 'k--', 'LineWidth', 3)
    plot(Dt*(1:N), xV_ekf(1,:), 'b-', 'LineWidth', 1)
    plot(Dt*(1:N), xV_ukf(1,:), 'r-', 'LineWidth', 1)
    plot(Dt*(1:N), xh(1,:), 'g-', 'LineWidth', 1)
    str = sprintf('Estimated state x_{1,k} vs. Time');
    title(str)
    ylabel('Pendulum angle (rad)')
    h_legend = legend('Meas', 'Real', 'EKF', 'UKF', 'PF');
    set(h_legend,'FontSize',10, 'Orientation', 'horizontal', 'Location', 'south');
    
    axis([0,9,-3,2])

%end
% subplot(2,1,2)
% plot(1:N, eV_ekf(1,:))
% title(sprintf('EKF RMSE vs. Time'))
% legend('EKF');
% axis([0,500,0,0.3])
%% Plot results
%figure
%for k=1:2                                 % plot results
%     subplot(2,3,3)
% %     figure
%      hold on
%     plot(1:N, sV_ukf(1,:), 'g-', 'LineWidth', 3)
%     plot(1:N, zV_ukf(1,:), 'r.', 'LineWidth', 1)
%     plot(1:N, xh(1,:), 'b-', 'LineWidth', 1)
% %     for i = 1:N
% %         hold on
% %         error_ellipse('C', blkdiag(PV_ukf{i}(1,1),1), 'mu', [i, xV_ukf(k,i)], 'style', 'r--')
% %         hold on
% %         error_ellipse('C', blkdiag(PV_ekf{i}(1,1),1), 'mu', [i, xV_ekf(k,i)], 'style', '--')
% %     end
%     str = sprintf('PF estimated state vs. Time');
%     title(str)
%     h_legend = legend('Real', 'Meas', 'PF');
%     set(h_legend,'FontSize',5);
%     axis([0,500,-4,4])

%end
subplot(2,3,[4:6])
hold on;
plot(Dt*(1:N), eV_obs(1,:), 'color', [0.7 0.7 0.7], 'LineWidth', 1);
plot(Dt*(1:N), eV_ekf(1,:), 'b-', 'LineWidth', 1)
plot(Dt*(1:N), eV_ukf(1,:), 'r-', 'LineWidth', 1)
plot(Dt*(1:N), err(1,:), 'g-', 'LineWidth', 1)
%plot(1:N, eV_ekf(1,:), 1:N, eV_ukf(1,:),1:N, err(1,:))
title(sprintf('RMSE vs. Time for x_{1,k}'))
h_legend = legend('Meas', 'EKF', 'UKF', 'PF');
set(h_legend,'FontSize',10);
xlabel('Time (s)')
ylabel('RMSE (rad)')
axis([0,9,0,0.25])