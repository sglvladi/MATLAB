% Number of Target Tracks
TrackNum = 2;

% Generate observations
lambdaV = 50;
lambda = lambdaV/(10^2); %1 0^2 is the area of the surveillance region
%[DataList,x1,y1] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, 0.1, 2, lambdaV, 1);
N=size(DataList,2) ; 

%% Initiate PF parameters
q=0.01;    %std of process 
r=0.1;    %std of measurement
nx=4;      %number of state dims

% Process equation x[k] = sys(k, x[k-1], u[k]);
sys_cch = @(k, xkm1, uk) [xkm1(1,:)+1*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+1*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)'];
obs = @(k, xk, vk) [xk(1)+vk(1); xk(2)+vk(2)];    % % Observation equation y[k] = obs(k, x[k], v[k]); (returns column vector)
nu = 4;                                           % size of the vector of process noise
gen_sys_noise_cch = @(u) mvnrnd(zeros(u, nu), diag([0,0,q^2,0.3^2])); 
nv = 2;                                           % size of the vector of observation noise
sigma_v = r;
cov_v = sigma_v^2*eye(nv);
p_obs_noise   = @(v) mvnpdf(v, zeros(1, nv), cov_v);
gen_obs_noise = @(v) mvnrnd(zeros(1, nv), cov_v);         % sample from p_obs_noise (returns column vector)
gen_x0_cch = @(Np) [10*rand(Np,1),10*rand(Np,1), mvnrnd(zeros(Np,1), 2*sigma_v^2), 2*pi*rand(Np,1)];
p_yk_given_xk = @(k, yk, xk) p_obs_noise((yk - obs(k, xk, zeros(1, nv)))');

% Assign PHD parameter values
par.k               = 1;                   % initial iteration number
par.Np              = 100000;              % number of particles
par.resampling_strategy = 'systematic_resampling'; % resampling strategy
par.sys = @(k, xkm1, uk) [xkm1(1,:)+1*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+1*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)']; % CH model
par.gen_x0 = @(Np) [10*rand(Np,1),10*rand(Np,1), mvnrnd(zeros(Np,1), 2*sigma_v^2), 2*pi*rand(Np,1)]; % Uniform position and heading, Gaussian speed
par.particles = par.gen_x0(par.Np)'; % Generate inital particles as per gen_x0
par.w = repmat(1/par.Np, par.Np, 1)'; % Uniform weights
par.likelihood = @(k, yk, xk) mvnpdf(yk, xk, cov_v); 
par.obs_model = @(xk) [xk(1,:); xk(2,:)];
par.clutter_flag = 1;
par.multi_flag = 1;
par.sys_noise = gen_sys_noise_cch;
par.Pbirth = 0.1;
par.Pdeath = 0.005;
par.J_k = 5000;
par.PD = 0.9;
par.lambda = lambda;
par.type = 'standard';

myphd = SMC_PHD(par);

figure
for i=1:N
    fprintf('Iteration = %d/\n',i);
    tempDataList = DataList{i}(:,:);
    tempDataList( :, ~any(tempDataList,1) ) = [];       
    
    % Change PHD filter parameters
    myphd.config.k = i; % Time index
    myphd.config.z = tempDataList; % New observations
    
    % Iterate PHD
    myphd.Predict();
    myphd.config.rhi = ones(1,size(myphd.config.z,2));
    myphd.Update();
    
    % Plot results
    clf;
    [bandwidth,density,X,Y]=kde2d(myphd.config.particles(1:2,:)');
    %contour3(X,Y,density,50);
    surf(X,Y,density);
    hold on;
    plot(myphd.config.particles(1,:), myphd.config.particles(2,:), '.')
    hold on;
    plot(myphd.config.z(1,:), myphd.config.z(2,:), 'y*');
    axis([0 10 0 10]);
    pause(0.01)
end