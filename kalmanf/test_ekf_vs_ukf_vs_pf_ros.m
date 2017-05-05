%% Clear everything
%clear
Dt = 1;
%% Initiate KF parameters
n=4;      %number of state
q=0.01;    %std of process 
r=0.3;    %std of measurement
s.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*10*q^2; % covariance of process
s.R=r^2*eye(n/2);        % covariance of measurement  
s.sys=@(t)(@(x)[x(1)+ t*x(3); x(2)+t*x(4); x(3); x(4)]);  % nonlinear state equations
s.obs=@(x)[x(1);x(2)];                               % measurement equation
st=[x_true(1);y_true(1)];                                % initial state
% Single-point initiation
Vmax = 0.4; % Max velocity = 0.4 m/s
s.x =[obs_x(1); obs_y(1); 0; 0]; %initial state
s.P = diag([q^2, q^2, (Vmax^2/3), (Vmax^2/3)]);                              % initial state covraiance
x_ukf = s.x;
P_ukf = s.P;
ukf_sys = @(x,u)[5*sin(x(2)+x(3))+u(1);2*cos(x(3))+u(2);x(3)+0.1+u(3)];
ukf_obs = @(x,v)[x(1);x(2);x(3)];

N=1019;%size(x,2);                                     % total dynamic steps

% Instantiate UKF and vars to store output
xV_ukf = zeros(n,N);          %estmate        % allocate memory
xV_ukf(:,1) = s.x;
PV_ukf = zeros(1,N);
%PV_ukf = cell(1,N);    % use to display ellipses 
sV_ukf = zeros(n/2,N);          %actual
sV_ukf(:,1) = st;
zV_ukf = zeros(n/2,N);
eV_ukf = [];
ukf = UKalmanFilter(s, 0.5, 0, 2);

% Instantiate EKF and vars to store output
xV_ekf = zeros(n,N);          %estmate        % allocate memory
xV_ekf(:,1) = s.x;
PV_ekf = zeros(1,N);    
%PV_ekf = cell(1,N);     % use to display ellipses 
sV_ekf = zeros(n/2,N);          %actual
sV_ekf(:,1) = st;
zV_ekf = zeros(n/2,N);
eV_ekf = [];
ekf = EKalmanFilter(s);

% Instantiate EKF and vars to store output
xV_kf = zeros(n,N);          %estmate        % allocate memory
xV_kf(:,1) = s.x;
PV_kf = zeros(1,N);    
%PV_ekf = cell(1,N);     % use to display ellipses 
sV_kf = zeros(n/2,N);          %actual
sV_kf(:,1) = st;
zV_kf = zeros(n/2,N);
eV_kf = [];
kf_s = s;
kf_s.sys = [1 0 Dt 0; 0 1 0 Dt; 0 0 1 0; 0 0 0 1];
kf_s.obs = [1 0 0 0; 0 1 0 0];
kf = KalmanFilter_new(kf_s);

%% Initiate PF parameters

% Process equation x[k] = sys(k, x[k-1], u[k]);
nx = 4;  % number of states
sys = @(k, xkm1, uk) [xkm1(1)+Dt*xkm1(3); xkm1(2)+Dt*xkm1(4); xkm1(3)+ uk(3); xkm1(4) + uk(4)]; % (returns column vector)

% Observation equation y[k] = obs(k, x[k], v[k]);
ny = 2;                                           % number of observations
obs = @(k, xk, vk) [xk(1)+vk(1); xk(2)+vk(2)];                  % (returns column vector)

% PDF of process noise and noise generator function
nu = 4;                                           % size of the vector of process noise
sigma_u = 0.01;
cov_u = [Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, 1]*10*sigma_u^2;
p_sys_noise   = @(u) mvnpdf(u, zeros(1, nu), cov_u);
gen_sys_noise = @(u) mvnrnd(zeros(1, nu), cov_u);         % sample from p_sys_noise (returns column vector)

% PDF of observation noise and noise generator function
nv = 2;                                           % size of the vector of observation noise
sigma_v = 0.3;
cov_v = sigma_v^2*eye(nv);
p_obs_noise   = @(v) mvnpdf(v, zeros(1, nv), cov_v);
gen_obs_noise = @(v) mvnrnd(zeros(1, nv), cov_v);         % sample from p_obs_noise (returns column vector)

% Initial PDF
% p_x0 = @(x) normpdf(x, 0,sqrt(10));             % initial pdf
gen_x0 = @(x) mvnrnd([obs_x(1); obs_y(1); 0; 0],diag([sigma_u^2, sigma_u^2, (Vmax^2/3), (Vmax^2/3)]));               % sample from p_x0 (returns column vector)

% Transition prior PDF p(x[k] | x[k-1])
% (under the suposition of additive process noise)
% p_xk_given_xkm1 = @(k, xk, xkm1) p_sys_noise(xk - sys(k, xkm1, 0));

% 
% (under the suposition of additive process noise)
p_yk_given_xk = @(k, yk, xk) p_obs_noise((yk - obs(k, xk, zeros(1, nv)))');

% Number of time steps
T = 1019;%1184;

% Separate memory space
%x = [x_true(:,1), y_true(:,1)]';
y = [obs_x'; obs_y']; % True state and observations
xh = zeros(nx, T); yh = zeros(ny, T); % Filtered state and observations

% Initial state and observatiopns
%xh0 = [obs_x(2); obs_y(2); obs_x(2)-obs_x(1); obs_y(2)-obs_y(1)] ; % 2-point state initialisation
% Single-point initiation
xh0=[obs_x(1); obs_y(1); 0; 0]; %initial state
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
pf.sys = sys;
pf.resampling_strategy = 'systematic_resampling';
my_pf = ParticleFilter(pf);


img = imread('maze.png');
 
% set the range of the axes
% The image will be stretched to this.
min_x = 0;
max_x = 10;
min_y = 0;
max_y = 10;
 
% make data to plot - just a line.
img_x = min_x:max_x;
img_y = (6/8)*x;

%% Track
figure
for k=2:N
    fprintf('Iteration = %d/%d\n',k,T);
    Dt = 1;
    % Compute process equation for KFs
    ukf.s.sys = s.sys(Dt);
    ekf.s.sys = s.sys(Dt);
    kf.s.sys = [1 0 Dt 0; 0 1 0 Dt; 0 0 1 0; 0 0 0 1];
    st = x(:,k);
    
    % Get next measurement
    ukf.s.z = [obs_x(k); obs_y(k)];                     % measurments
    ekf.s.z = ukf.s.z;
    kf.s.z = ukf.s.z;
    my_pf.pf.k = k;
    my_pf.pf.z = y(:,k);

    % Store new state and measurement
    sV_ukf(:,k)= x(:,k);                             % save actual state
    zV_ukf(:,k)  = ukf.s.z;                             % save measurment
    sV_ekf(:,k)= st;                             % save actual state
    zV_ekf(:,k)  = ekf.s.z; 

    % Iterate both filters
    kf.s = kf.Iterate(kf.s);
    ekf.s = ekf.Iterate(ekf.s);
    ukf.s = ukf.Iterate(ukf.s);            % ekf 
    my_pf.pf = my_pf.Iterate(my_pf.pf);


    %% Store estimated state and covariance
    xV_ukf(:,k) = ukf.s.x;%x_ukf;                            % save estimate
    PV_ukf(k)= ukf.s.P(1,1); % P_ukf(1,1) ; 
    %PV_ukf{k}= ukf.s.P;    % Use to store whole covariance matrix
    xV_ekf(:,k) = ekf.s.x;                            % save estimate
    PV_ekf(k) = ekf.s.P(1,1);
    %PV_ekf{k} = ekf.s.P;    % Use to store whole covariance matrix
    xV_kf(:,k) = kf.s.x;                            % save estimate
    PV_kf(k) = kf.s.P(1,1); 
    
    xh(:,k) = my_pf.pf.xhk(:,k);
    % filtered observation
    yh(:,k) = obs(k, xh(:,k), zeros(1, nv));

    %% Compute squared error
%     eV_ukf(:,k) = sum((ukf.s.x(:,1) - [st;v(:,k)]).*(ukf.s.x(:,1) - [st;v(:,k)]),1)/4
%     eV_ekf(:,k) = sum((ekf.s.x(:,1) - [st;v(:,k)]).*(ekf.s.x(:,1) - [st;v(:,k)]),1)/4
%     eV_kf(:,k) = sum((kf.s.x(:,1) - [st;v(:,k)]).*(kf.s.x(:,1) - [st;v(:,k)]),1)/4

    clf;
    imagesc([min_x max_x], [min_y max_y], flipud(img));
    hold on;
    h1 = plot(sV_ukf(1,1:k),sV_ukf(2,1:k),'k--','LineWidth',1);
    h2 = plot(sV_ukf(1,k),sV_ukf(2,k),'ko','MarkerSize', 20);
    h3 = plot(zV_ukf(1,1:k), zV_ukf(2,1:k),'k.','LineWidth',1);
    h4 = plot(xV_kf(1,1:k),xV_kf(2,1:k),'c','LineWidth',1);
    h5 = plot(xV_kf(1,k),xV_kf(2,k),'co','MarkerSize', 20);
    h6 = plot(xV_ekf(1,1:k),xV_ekf(2,1:k),'b','LineWidth',1);
    h7 = plot(xV_ekf(1,k),xV_ekf(2,k),'bo','MarkerSize', 20);
    h8 = plot(xV_ukf(1,1:k),xV_ukf(2,1:k),'r','LineWidth',1);
    h9 = plot(xV_ukf(1,k),xV_ukf(2,k),'ro','MarkerSize', 20);
    h10 = plot(xh(1,1:k),xh(2,1:k),'g','LineWidth',1);
    h11 = plot(xh(1,k),xh(2,k),'go','MarkerSize', 20);
    legend([h1 h3 h4 h6 h8 h10],'Ground Truth', 'measurements', 'KF', 'EKF', 'UKF', 'PF');
    title('State vs estimated state by the particle filter vs particle paths','FontSize',14);
    %set the y-axis back to normal.
    set(gca,'ydir','normal');
    pause(0.001)
    % Generate new state
   % st = ukf.s.sys(st)+q*(-1 + 2*rand(3,1));                % update process 
end

%% Compute & Print RMSE
eV_kf = (xV_kf(1:2,:) - x(:,1:N)).*(xV_kf(1:2,:) - x(:,1:N));
eV_ekf = (xV_ekf(1:2,:) - x(:,1:N)).*(xV_ekf(1:2,:) - x(:,1:N));
eV_ukf = (xV_ukf(1:2,:) - x(:,1:N)).*(xV_ukf(1:2,:) - x(:,1:N));
RMSE_kf = sqrt(sum(eV_kf,2)/N)
RMSE_ekf = sqrt(sum(eV_ekf,2)/N)
RMSE_ukf = sqrt(sum(eV_ukf,2)/N)

%% Compute velocity
vV_ukf = [];
vV_ukf(:,1) = [zV_ukf(1,2)-zV_ukf(1,1); zV_ukf(2,2)-zV_ukf(2,1)];
for i=2:size(zV_ukf, 2)
    vV_ukf(:,i) = [zV_ukf(1,i)-zV_ukf(1,i-1); zV_ukf(2,i)-zV_ukf(2,i-1)];
end
%% Compute RMSE
pf_err = (xh(1:2,:) - x(:,1:N)).*(xh(1:2,:) - x(:,1:N));
err = (zV_ukf(:,:) - x(:,1:N)).*(zV_ukf(:,:) - x(:,1:N));
RMSE_pf = sqrt(sum(pf_err,2)/T)
RMSE_meas = sqrt(sum(err,2)/T)
%% Plot results
figure('units','centimeters','position',[.1 .1 10.2 9.1])
for k=1:4                                 % plot results
%     figure
    subplot(4,1,k)
    %figure
    if (k<3)
        axis([0,100,0,0.18])
    elseif (k==3)
        axis([0,100,0,0.06])
    else
        axis([0,100,0,0.1])
    end
     hold on
    plot(1:N, sqrt(err(k,:)), 'k--', 'LineWidth', 1) 
    plot(1:N, sqrt(eV_kf(k,:)), 'c-', 'LineWidth', 1)
    plot(1:N, sqrt(eV_ekf(k,:)), 'b-', 'LineWidth', 1)
    plot(1:N, sqrt(eV_ukf(k,:)), 'r-', 'LineWidth', 1)
    plot(1:N, sqrt(pf_err(k,:)), 'g-', 'LineWidth', 1) 
    %xlabel('label_{subscript}')
    str = sprintf(['RMSE vs. Time for x_{' num2str(k) ',k}']);
    title(str)
    h_legend = legend('Meas', 'KF', 'EKF', 'UKF', 'PF','Orientation','horizontal');
    set(h_legend,'FontSize',9);
    %text(0,0.90,'$\textcircled{a}$', 'Interpreter', 'latex');


end
%% Plot results
% figure
% for k=1:2                                 % plot results
%     subplot(3,1,k)
% %     figure
% %     hold on
%     plot(1:N, sV_ukf(k,:), 'k--', 1:N, xV_ukf(k,:), 'b-',1:N, xV_ekf(k,:), 'g-', 1:N, zV_ukf(k,:), 'r.')
% %     for i = 1:N
% %         hold on
% %         error_ellipse('C', blkdiag(PV_ukf{i}(1,1),1), 'mu', [i, xV_ukf(k,i)], 'style', 'r--')
% %         hold on
% %         error_ellipse('C', blkdiag(PV_ekf{i}(1,1),1), 'mu', [i, xV_ekf(k,i)], 'style', '--')
% %     end
%     str = sprintf('EKF vs UKF estimated state X(%d)',k);
%     title(str)
%     legend('Real', 'UKF', 'EKF', 'Meas');
% end
% subplot(3,1,3)
% plot(1:N, PV_ukf(:,:), 1:N, PV_ekf(:,:))
% title(sprintf('EKF vs UKF estimated covariance P(1,1)',k))
% legend('UKF', 'EKF');
% 
% figure
% plot(sV_ukf(2,:), eV_ukf(:),xV_ukf(1,:), xV_ukf(2,:), xV_ekf(1,:), xV_ekf(2,:));
% title(sprintf('EKF vs UKF estimated covariance P(1,1)',k))
% legend('true', 'UKF', 'EKF');
% 
% figure
% plot(sV_ukf(1,:), sV_ukf(2,:),xV_ukf(1,:), xV_ukf(2,:), xV_ekf(1,:), xV_ekf(2,:));
% title(sprintf('EKF vs UKF estimated covariance P(1,1)',k))
% legend('true', 'UKF', 'EKF');

