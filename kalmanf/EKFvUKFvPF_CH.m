PV_err = []
No_of_Runs = 10
r_list = [1]
for r_iter = 1:length(r_list)
    for Run_iter=1:No_of_Runs
        [x, obs_x, obs_y, obs_r, obs_phi] = gen_obs(x_true, y_true, th_true, r_list(r_iter));
        %% Clear everything
        clear pf;
        clear F;
        Dt = 1;
        %% Global parameters
        nx=4;      %number of state dims
        ny = 2;    %number of observation dims
        q=0.01;    %std of process 
        % 120 210-230 270 346-360 -370-390-460
        r=r_list(r_iter);    %std of measurement                               % measurement equation
        st=[x_true(1);y_true(1)];                                % initial state
        % Number of iterations
        N=size(x,2);                                     % total dynamic steps
        T=N;

        % Initial state and observatiopns
        %s.x = [obs_x(2); obs_y(2); obs_x(2)-obs_x(1); obs_y(2)-obs_y(1)] ; % 2-point state initialisation
        % Single-point initiation
        Vmax = 0.4; % Max velocity = 0.4 m/s

        %% KF bank
        s.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2; % covariance of process
        s.R=r^2*eye(nx/2);        % covariance of measurement  
        s.sys=@(t)(@(x)[x(1)+ t*x(3); x(2)+t*x(4); x(3); x(4)]);  % nonlinear state equations
        s.obs=@(x)[x(1);x(2)];
        s.obs2=@(x,u)[x(1)+u(1);x(2)+u(2)];
        s.x =[x_true(1); y_true(1); x_true(2)-x_true(1); y_true(2)-y_true(1)]; %initial state
        s.P = diag([q^2, q^2, (Vmax^2/3), 0.3^2]);                              % initial state covraiance

        %% Estimated State container
        xh = zeros(6,nx, T); xh(1,:,1) = s.x; xh(2,:,1) = s.x; xh(3,:,1) = s.x; xh(4,:,1) = s.x; xh(5,:,1) = s.x; xh(6,:,1) = s.x; xh(7,:,1) = s.x;
        %yh = zeros(5,ny, T); yh(1,:,1) = obs(1, xh0, zeros(1, nv)); yh(2,:,1) = obs(1, xh0, zeros(1, nv)); yh(3,:,1) = obs(1, xh0, zeros(1, nv)); yh(4,:,1) = obs(1, xh0, zeros(1, nv));
        pV_err = [];
        hV_err = [];

        %% Instantiate UKF and vars to store output
        %s.sys=@(t)(@(x,u)[x(1)+ t*x(3)*cos(x(4))+u(1); x(2)+t*x(3)*sin(x(4))+u(2); x(3)+u(3); x(4)+u(4)]);
        s.sys=@(t)(@(x,u)[x(1)+ t*x(3)+u(1); x(2)+t*x(4)+u(2); x(3)+u(3); x(4)+u(4)]);
        %s.sys=s.sys3;
        s.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2;
        %s.Q = diag([0.01^10,0.01^10,0.01^2,0.3^2]);
        s.obs = s.obs2;
        sV_ukf = zeros(nx/2,N);          %actual
        sV_ukf(:,1) = st;
        zV_ukf = zeros(nx/2,N);
        % eV_ukf = [];
        ukf_ccv = UKalmanFilter(s, 0.5, 0, 2);

        s.sys=@(t)(@(x,u)[x(1)+ t*x(3)*cos(x(4))+u(1); x(2)+t*x(3)*sin(x(4))+u(2); x(3)+u(3); x(4)+u(4)]);
        s.x =[x_true(1); y_true(1); sqrt((x_true(2)-x_true(1))^2+(y_true(2)-y_true(1))^2); atan((y_true(2)+0.1-y_true(1))/(x_true(2)+0.1-x_true(1)))];
        s.P = diag([q^2, q^2, (Vmax^2/3), 0.3^2]);
        %s.sys=@(t)(@(x,u)[x(1)+ t*x(3)+u(1); x(2)+t*x(4)+u(2); x(3)+u(3); x(4)+u(4)]);
        s.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2;
        s.sys3=s.sys;
        s.Q = diag([0.01^10,0.01^10,0.01^2,0.3^2]);
        s.obs = s.obs2;
        % eV_ukf = [];
        ukf = UKalmanFilter(s, 0.5, 0, 2);

        %% Instantiate EKF and vars to store output
        s.sys=@(t)(@(x)[x(1)+ t*x(3)*cos(x(4)); x(2)+t*x(3)*sin(x(4)); x(3); x(4)]);  % nonlinear state equations
    %     s.sys=@(t)(@(x)[x(1)+ t*x(3); x(2)+t*x(4); x(3); x(4)]);
        s.obs=@(x)[x(1);x(2)];
        s.Q = diag([0,0,0.01^2,0.3^2]);
        %s.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2;
        s.P = diag([q^2, q^2, (Vmax^2/3), 0.3^2]);
        ekf = EKalmanFilter(s);
        %ukf = UKalmanFilter(s, 0.5, 0, 2);
        s.sys2 = s.sys;
        s.Q2 = s.Q;
        s.sys=@(t)(@(x)[x(1)+ t*x(3); x(2)+t*x(4); x(3); x(4)]);  % nonlinear state equations
        s.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2; % covariance of process

        %% Instantiate KF and vars to store output
        s.x =[x_true(1); y_true(1); x_true(2)-x_true(1); y_true(2)-y_true(1)];
        s.P = diag([q^2, q^2, (Vmax^2/3), (Vmax^2/3)]);
        kf_s = s;
        kf_s.sys = [1 0 Dt 0; 0 1 0 Dt; 0 0 1 0; 0 0 0 1];
        kf_s.obs = [1 0 0 0; 0 1 0 0];
        kf = KalmanFilter_new(kf_s);

        %% Initiate PF parameters

        % Process equation x[k] = sys(k, x[k-1], u[k]);
        sys_pchr = @(k, xkm1, uk) [sqrt((xkm1(1)^2 +2*xkm1(1)*xkm1(3)*Dt*cos(xkm1(2)-xkm1(4))+xkm1(3)^2*Dt^2)); atan((xkm1(1)*sin(xkm1(2))+xkm1(3)*Dt*sin(xkm1(4)))/(xkm1(1)*cos(xkm1(2))+xkm1(3)*Dt*cos(xkm1(4)))); xkm1(3)+ uk(3); xkm1(4) + xkm1(5); xkm1(5) + uk(5)]; % (returns column vector)
        sys_pchr2 = @(k, xkm1, uk) [sqrt((xkm1(1)*cos(xkm1(2))+xkm1(3)*Dt*cos(xkm1(4)))^2+(xkm1(1)*sin(xkm1(2))+xkm1(3)*Dt*sin(xkm1(4)))^2); atan((xkm1(1)*sin(xkm1(2))+xkm1(3)*Dt*sin(xkm1(4)))/(xkm1(1)*cos(xkm1(2))+xkm1(3)*Dt*cos(xkm1(4)))); xkm1(3)+ uk(3); xkm1(4) + xkm1(5); xkm1(5) + uk(5)]; % (returns column vector)
        sys_pch = @(k, xkm1, uk) [sqrt((xkm1(1)*cos(xkm1(2))+xkm1(3)*Dt*cos(xkm1(4)))^2+(xkm1(1)*sin(xkm1(2))+xkm1(3)*Dt*sin(xkm1(4)))^2); atan((xkm1(1)*sin(xkm1(2))+xkm1(3)*Dt*sin(xkm1(4)))/(xkm1(1)*cos(xkm1(2))+xkm1(3)*Dt*cos(xkm1(4)))); xkm1(3)+ uk(3); xkm1(4) + uk(4)]; % (returns column vector)
        sys_ccv = @(k, xkm1, uk) [xkm1(1,:)+1*xkm1(3,:)+uk(:,1)'; xkm1(2,:)+1*xkm1(4,:)+uk(:,2)'; xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)'];
        sys_cch = @(k, xkm1, uk) [xkm1(1,:)+1*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+1*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)'];
        sys_cchr = @(k, xkm1, uk) [xkm1(1)+1*xkm1(3)*cos(xkm1(4)); xkm1(2)+1*xkm1(3)*sin(xkm1(4)); xkm1(3)+ uk(3); xkm1(4) + xkm1(5); xkm1(5) + uk(5)]; % (returns column vector)

        % Observation equation y[k] = obs(k, x[k], v[k]);
        obs = @(k, xk, vk) [xk(1,:)+vk(1,:); xk(2,:)+vk(2,:)];                  % (returns column vector)
        obs_p = @(k, xk, vk) [xk(1)*cos(xk(2))+vk(1); xk(1)*sin(xk(2))+vk(2)];
        obs_radar = @(k, xk, vk) [sqrt(xk(1)^2+xk(2)^2); atan2(xk(2),xk(1))];

        % PDF of process noise and noise generator function
        nu = 4;                                           % size of the vector of process noise
        sigma_u = q;
        cov_u = [Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2;
        gen_sys_noise_pch = @(u) mvnrnd(zeros(size(u,2), nu), diag([0,0,0.01^2,0.3^2]));
        gen_sys_noise_pchr = @(u) mvnrnd(zeros(size(u,2), 5), diag([0,0, 0.01^2, 0, 0.1^2]));         % sample from p_sys_noise (returns column vector)
        gen_sys_noise_ccv = @(u) mvnrnd(zeros(size(u,2), nu), 5*cov_u);
        gen_sys_noise_ccv2 = @(u) mvnrnd(zeros(size(u,2), nu), 2*cov_u);
        gen_sys_noise_cch = @(u) mvnrnd(zeros(size(u,2), nu), diag([0,0,0.01^2,0.3^2])); 
        gen_sys_noise_cchr = @(u) mvnrnd(zeros(size(u,2),5), diag([0,0,0.01^2,0,1^2]));
        % PDF of observation noise and noise generator function
        nv = 2;                                           % size of the vector of observation noise
        sigma_v = r;
        cov_v = sigma_v^2*eye(nv);
        p_obs_noise   = @(v) mvnpdf(v, zeros(1, nv), cov_v);
        gen_obs_noise = @(v) mvnrnd(zeros(size(v,2), nv), cov_v);         % sample from p_obs_noise (returns column vector)
        p_obs_noise_radar   = @(v) mvnpdf(v, zeros(1, nv), [r,0;0,2*pi/360]);
        gen_obs_noise_radar = @(v) mvnrnd(zeros(1, nv), [r,0;0,2*pi/360]);         % sample from p_obs_noise (returns column vector)

        % Initial PDF
        gen_x0_ccv = @(Np) mvnrnd(repmat([x_true(1), y_true(1), x_true(2)-x_true(1),y_true(2)-y_true(1)], Np, 1),diag([sigma_u^2, sigma_u^2, (Vmax^2/3), (Vmax^2/3)]));               % sample from p_x0 (returns column vector)              % sample from p_x0 (returns column vector)
        gen_x0_cch = @(Np) mvnrnd(repmat([x_true(1), y_true(1), sqrt((x_true(2)-x_true(1))^2+(y_true(2)-y_true(1))^2), atan((y_true(2)+0.1-y_true(1))/(x_true(2)+0.1-x_true(1)))], Np, 1),diag([sigma_u^2, sigma_u^2, (Vmax^2/3), 0.3^2]));
        gen_x0_cchr = @(Np) mvnrnd(repmat([obs_x(1), obs_y(1), 0, 0, 0], Np , 1),diag([sigma_u^2, sigma_u^2, (Vmax^2/3), 0, 0]));
        gen_x0_pch = @(Np) mvnrnd(repmat([sqrt(obs_x(1)^2+obs_y(1)^2), atan(obs_y(1)/obs_x(1)), 0, 0], Np, 1),diag([sigma_u^2, sigma_u^2, (Vmax^2/3), 0]));
        gen_x0_pchr = @(Np) mvnrnd(repmat([sqrt(obs_x(1)^2+obs_y(1)^2), atan(obs_y(1)/obs_x(1)), 0, 0, 0], Np, 1),diag([sigma_u^2, sigma_u^2, (Vmax^2/3), 0, 0]));

        % Observation likelihood PDF p(y[k] | x[k])
        % (under the suposition of additive process noise)
        p_yk_given_xk = @(k, yk, xk) p_obs_noise((yk*ones(1, size(xk,2)) - obs(k, xk, zeros(nv,size(xk,2))))');
        p_yk_given_xk_p = @(k, yk, xk) p_obs_noise((yk - obs_p(k, xk, zeros(1, nv)))');
        p_yk_given_xk_radar = @(k, yk, xk) p_obs_noise_radar((yk - obs_radar(k, xk, zeros(1, nv)))');

        % Separate memory space
        %x = [x_true(:,1), y_true(:,1)]';
        y = [obs_x'; obs_y']; % True state and observations
        y_radar = [obs_r'; obs_phi']; % True state and observations
        % Assign PF parameter values
        pf.clutter_flag = 0;
        pf.k               = 1;                   % initial iteration number
        pf.Np              = 10000;                 % number of particles
        %pf.w               = zeros(pf.Np, T);     % weights
        pf.particles       = zeros(5, pf.Np); % particles
        pf.gen_x0          = gen_x0_pchr;              % function for sampling from initial pdf p_x0 (used to generate initial particles)
        pf.obs             = p_yk_given_xk_p;       % function of the observation likelihood PDF p(y[k] | x[k])
        pf.sys_noise       = gen_sys_noise_pchr;       % function for generating system noise
        %pf.p_x0 = p_x0;                          % initial prior PDF p(x[0])
        %pf.p_xk_given_ xkm1 = p_xk_given_xkm1;   % transition prior PDF p(x[k] | x[k-1])
        pf.xhk = s.x;
        pf.sys = sys_pchr;
        pf.resampling_strategy = 'systematic_resampling';

        % PF-PCHR
        %pf_pchr = ParticleFilterMin(pf);

        % PF-CCV
        pf.sys = sys_ccv;
        pf.particles = zeros(nx, pf.Np); % particles
        pf.gen_x0 = gen_x0_ccv;
        pf.obs = p_yk_given_xk;
        pf.sys_noise = gen_sys_noise_ccv;
        %pf.kf = ukf;
        %pf.kf.s.sys = s.sys3(Dt);
        pf_ccv = ParticleFilterMin2(pf);

        % PF-CCH
        pf.sys = sys_cch;
        pf.particles = zeros(nx, pf.Np); % particles
        pf.gen_x0 = gen_x0_cch;
        pf.obs = p_yk_given_xk;
        pf.sys_noise = gen_sys_noise_cch;
        %pf.kf = ukf;
        %pf.kf.s.sys = s.sys3(Dt);
        pf_cch = ParticleFilterMin2(pf);

        % PF-CCV-OPT
        pf.sys = sys_ccv;
        pf.particles = zeros(nx, pf.Np); % particles
        pf.gen_x0 = gen_x0_ccv;
        pf.obs = p_yk_given_xk;
        pf.sys_noise = gen_sys_noise_ccv;
        pf.kf = ukf_ccv;
        pf.kf.s.sys = pf.kf.s.sys(Dt);
        pf_ccv_opt = ParticleFilterMin2(pf);

        % PF-CCH-OPT
        pf.sys = sys_cch;
        pf.particles = zeros(nx, pf.Np); % particles
        pf.gen_x0 = gen_x0_cch;
        pf.obs = p_yk_given_xk;
        pf.kf = ukf;
        pf.kf.s.sys = pf.kf.s.sys(Dt);
        pf.sys_noise = gen_sys_noise_cch;
        pf_cch_opt = ParticleFilterMin2(pf);

        % PF-PCH
    %     pf.sys = sys_pch;
    %     pf.particles       = zeros(nx, pf.Np); % particles
    %     pf.gen_x0 = gen_x0_pch;
    %     pf.obs             = p_yk_given_xk_p; 
    %     pf.sys_noise = gen_sys_noise_pch;
    %     pf_pch = ParticleFilterMin(pf);


        % img = imread('maze.png');
        %  
        % % set the range of the axes
        % % The image will be stretched to this.
        % min_x = -5;
        % max_x = 15;
        % min_y = -5;
        % max_y = 15;
        %  
        % % make data to plot - just a line.
        % img_x = min_x:max_x;
        % img_y = (6/8)*x;

        %% Track
        %figure
        %figure
        for k=2:N
            fprintf('R = %d/%d\n',r_iter,length(r_list));
            fprintf('Run = %d/%d\n',Run_iter,No_of_Runs);
            fprintf('Iteration = %d/%d\n',k,T);
            Dt = 1;
            % Compute process equation for KFs
            ukf.s.sys = s.sys3(Dt);
            ekf.s.sys = s.sys2(Dt);
            kf.s.sys = [1 0 Dt 0; 0 1 0 Dt; 0 0 1 0; 0 0 0 1];
            st = x(1:2,k);

            % Get next measurement
            ukf.s.z = [obs_x(k); obs_y(k)];                     % measurments
            ekf.s.z = ukf.s.z;
            kf.s.z = ukf.s.z;
            %pf_pchr.pf.k = k;
            %pf_pchr.pf.z = y(:,k);
            pf_ccv.pf.k = k;
            pf_ccv.pf.z = ukf.s.z;
            pf_cch.pf.k = k;
            pf_cch.pf.z = ukf.s.z;
            pf_ccv_opt.pf.k = k;
            pf_ccv_opt.pf.z = ukf.s.z;
            pf_cch_opt.pf.k = k;
            pf_cch_opt.pf.z = ukf.s.z;

            % Store new state and measurement
            sV_ukf(:,k)= x(1:2,k);                             % save actual state
            zV_ukf(:,k)  = ukf.s.z;                             % save measurment
        %     sV_ekf(:,k)= st;                             % save actual state
        %     zV_ekf(:,k)  = ekf.s.z; 

            % Iterate both filters
            kf.s = kf.Iterate(kf.s);
            ekf.s = ekf.Iterate(ekf.s);
            ukf.s = ukf.Iterate(ukf.s);            % ekf 
            %tic
            %pf_pchr.pf = pf_pchr.Iterate(pf_pchr.pf);
            %toc
            %tic
            pf_ccv.pf = pf_ccv.Iterate(pf_ccv.pf);
            pf_cch.pf = pf_cch.Iterate(pf_cch.pf);
            %toc
            %tic
            pf_ccv_opt.pf = pf_ccv_opt.Iterate(pf_ccv_opt.pf);
            %toc
            %tic
            pf_cch_opt.pf = pf_cch_opt.Iterate(pf_cch_opt.pf);
            %toc

            %% Store estimated state and covariance
        %     xV_ukf(:,k) = ukf.s.x;%x_ukf;                            % save estimate
        %     PV_ukf(k)= ukf.s.P(1,1); % P_ukf(1,1) ; 
        %     %PV_ukf{k}= ukf.s.P;    % Use to store whole covariance matrix
        %     xV_ekf(:,k) = ekf.s.x;                            % save estimate
        %     PV_ekf(k) = ekf.s.P(1,1);
        %     %PV_ekf{k} = ekf.s.P;    % Use to store whole covariance matrix
        %     xV_kf(:,k) = kf.s.x;                            % save estimate
        %     PV_kf(k) = kf.s.P(1,1); 

            xh(1,:,k) = kf.s.x;
            pV_err(1,k) = sqrt((xh(1,1,k) - x(1,k))^2 + (xh(1,2,k)-x(2,k))^2);
            hV_err(1,k) = (atan(xh(1,3,k)/xh(1,4,k))-x(3,k)); 
            xh(2,:,k) = ekf.s.x;
            pV_err(2,k) = sqrt((xh(2,1,k) - x(1,k))^2 + (xh(2,2,k)-x(2,k))^2);
            hV_err(2,k) = (rem(xh(2,4,k),pi)-x(3,k));
            xh(3,:,k) = ukf.s.x;
            pV_err(3,k) = sqrt((xh(3,1,k) - x(1,k))^2 + (xh(3,2,k)-x(2,k))^2);
            hV_err(3,k) = (atan(xh(3,3,k)/xh(3,4,k))-x(3,k)); 
            xh(4,:,k) = pf_ccv.pf.xhk(:);
            pV_err(4,k) = sqrt((xh(4,1,k) - x(1,k))^2 + (xh(4,2,k)-x(2,k))^2);
            hV_err(4,k) = (rem(xh(4,4,k),pi)-x(3,k));
            xh(5,:,k) = pf_cch.pf.xhk(:);
            pV_err(5,k) = sqrt((xh(5,1,k) - x(1,k))^2 + (xh(5,2,k)-x(2,k))^2);
            hV_err(5,k) = (rem(xh(5,4,k),pi)-x(3,k));
            xh(6,:,k) = pf_ccv_opt.pf.xhk(1:4);
            pV_err(6,k) = sqrt((xh(6,1,k) - x(1,k))^2 + (xh(6,2,k)-x(2,k))^2);
            hV_err(6,k) = (rem(xh(6,4,k),pi)-x(3,k));
            xh(7,:,k) = pf_cch_opt.pf.xhk(:);
            pV_err(7,k) = sqrt((xh(7,1,k) - x(1,k))^2 + (xh(7,2,k)-x(2,k))^2);
            hV_err(7,k) = (rem(xh(7,4,k),pi)-x(3,k));
            for i=1:size(hV_err,1)
                if hV_err(i,k)<0
                    hV_err(i,k) = hV_err(i,k) + pi;
                end
            end
            pV_err(8,k) = sqrt((obs_x(k,1) - x(1,k))^2 + (obs_y(k,1)-x(2,k))^2);
        %     xh(:,k) = my_pf.pf.xhk(:,k);
        %     % filtered observation
        %     yh(:,k) = obs(k, xh(:,k), zeros(1, nv));

            %% Compute squared error
        %     eV_ukf(:,k) = sum((ukf.s.x(:,1) - [st;v(:,k)]).*(ukf.s.x(:,1) - [st;v(:,k)]),1)/4
        %     eV_ekf(:,k) = sum((ekf.s.x(:,1) - [st;v(:,k)]).*(ekf.s.x(:,1) - [st;v(:,k)]),1)/4
        %     eV_kf(:,k) = sum((kf.s.x(:,1) - [st;v(:,k)]).*(kf.s.x(:,1) - [st;v(:,k)]),1)/4
        %std(pf_cch.pf.particles(1,:))
        %std(pf_cch.pf.particles(1,:),pf_cch.pf.w')
    %         clf;
    %         axis([-5 30 -5 30])
    %     %     imagesc([min_x max_x], [min_y max_y], flipud(img));
    %          hold on;
    %          h1 = plot(sV_ukf(1,1:k),sV_ukf(2,1:k),'k--','LineWidth',1);
    %          h2 = plot(sV_ukf(1,k),sV_ukf(2,k),'ko','MarkerSize', 20);
    %          h3 = plot(zV_ukf(1,1:k), zV_ukf(2,1:k),'k.','LineWidth',1);
    % %         h3 = plot(obs_r(1:k).*cos(obs_phi(1:k)), obs_r(1:k).*sin(obs_phi(1:k)),'r.','LineWidth',1);
    %         h4 = plot(permute(xh(1,1,1:k),[2 3 1]),permute(xh(1,2,1:k),[2 3 1]),'c','LineWidth',1);
    %         h5 = plot(permute(xh(1,1,k),[2 3 1]),permute(xh(1,2,k),[2 3 1]),'co','MarkerSize', 20);
    %          h6 = plot(permute(xh(2,1,1:k),[2 3 1]),permute(xh(2,2,1:k),[2 3 1]),'b','LineWidth',1);
    %          h7 = plot(permute(xh(2,1,k),[2 3 1]),permute(xh(2,2,k),[2 3 1]),'bo','MarkerSize', 20);
    %         % h2=plot_gaussian_ellipsoid(ekf.s.x(1:2,1), ekf.s.P(1:2,1:2));
    %         h8 = plot(permute(xh(3,1,1:k),[2 3 1]),permute(xh(3,2,1:k),[2 3 1]),'r','LineWidth',1);
    %         h9 = plot(permute(xh(3,1,k),[2 3 1]),permute(xh(3,2,k),[2 3 1]),'ro','MarkerSize', 20);
    %         h10 = plot(permute(xh(4,1,1:k),[2 3 1]),permute(xh(4,2,1:k),[2 3 1]),'Color',[1,0.6,0],'LineWidth',1);
    %         h11 = plot(permute(xh(4,1,k),[2 3 1]),permute(xh(4,2,k),[2 3 1]),'Marker','o','Color',[1,0.6,0],'MarkerSize', 20);
    %         %h2=plot_gaussian_ellipsoid([permute(xh(4,1,k),[2 3 1]);permute(xh(4,2,k),[2 3 1])], [std(pf_ccv.pf.particles(1,:))^2,0;0,std(pf_ccv.pf.particles(2,:))^2]);
    %          h12 = plot(permute(xh(5,1,1:k),[2 3 1]),permute(xh(5,2,1:k),[2 3 1]),'m','LineWidth',1);
    %         %plot(pf_ccv_opt.pf.particles(1,:),pf_ccv_opt.pf.particles(2,:),'.')
    %          %fit_ellipse(pf_cch.pf.particles(1,:),pf_cch.pf.particles(2,:), gca)
    %         % h2=plot_gaussian_ellipsoid([permute(xh(5,1,k),[2 3 1]);permute(xh(5,2,k),[2 3 1])], [std(pf_cch.pf.particles(1,:))^2,0;0,std(pf_cch.pf.particles(2,:))^2]);
    %          h13 = plot(permute(xh(5,1,k),[2 3 1]),permute(xh(5,2,k),[2 3 1]),'mo','MarkerSize', 20);
    %         h14 = plot(permute(xh(6,1,1:k),[2 3 1]),permute(xh(6,2,1:k),[2 3 1]),'g','LineWidth',1);
    %         h15 = plot(permute(xh(6,1,k),[2 3 1]),permute(xh(6,2,k),[2 3 1]),'go','MarkerSize', 20);
    %         %plot(pf_ccv.pf.particles(1,:),pf_ccv.pf.particles(2,:), 'r.');
    %         %h2=plot_gaussian_ellipsoid([permute(xh(6,1,k),[2 3 1]);permute(xh(6,2,k),[2 3 1])], [std(pf_ccv_opt.pf.particles(1,:))^2,0;0,std(pf_ccv_opt.pf.particles(2,:))^2]);
    %         h16 = plot(permute(xh(7,1,1:k),[2 3 1]),permute(xh(7,2,1:k),[2 3 1]),'Color',[0.8,0.8,0],'LineWidth',1);
    %         h17 = plot(permute(xh(7,1,k),[2 3 1]),permute(xh(7,2,k),[2 3 1]),'Marker','o','Color',[0.8,0.8,0],'MarkerSize', 20);
    %         %h2=plot_gaussian_ellipsoid([permute(xh(7,1,k),[2 3 1]);permute(xh(7,2,k),[2 3 1])], [std(pf_cch_opt.pf.particles(1,:))^2,0;0,std(pf_cch_opt.pf.particles(2,:))^2]);
    %         % %         h16 = plot(permute(xh(7,1,1:k),[2 3 1]),permute(xh(7,2,1:k),[2 3 1]),'Color',[1,0.6,0.4],'LineWidth',1);
    % % %         h17 = plot(permute(xh(7,1,k),[2 3 1]),permute(xh(7,2,k),[2 3 1]),'Marker','o','Color',[1,0.6,0.4],'MarkerSize', 20);
    %          legend([h1 h3 h4 h6 h8 h10 h12 h14 h16],'Ground Truth', 'measurements', 'KF-CCV', 'EKF-CCH', 'UKF-CCH', 'PF-CCV', 'PF-CCH', 'PF-CCV-OPT', 'PF-CCH-OPT');%, 'PF-PCH', 'PF-PCHR');
    %          title('State vs estimated state by the particle filter vs particle paths','FontSize',14);
    % %         %set the y-axis back to normal.
    % %         set(gca,'ydir','normal');
    %          pause(0.001)
            % Generate new state
    %         ax = gca;
    %         ax.Units = 'pixels';
    %         pos = ax.Position;
    %         marg = 30;
    %         rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
    %         F(k) = getframe(gcf);
    %         ax.Units = 'normalized';
           % st = ukf.s.sys(st)+q*(-1 + 2*rand(3,1));                % update process 
        end
    %     F = F(2:end);
    %     name = sprintf('test_s%d.avi',Run_iter);
    %     vidObj = VideoWriter(name);
    %     vidObj.Quality = 100;
    %     vidObj.FrameRate = 10;
    %     open(vidObj);
    %     writeVideo(vidObj, F);
    %     close(vidObj);
    %     winopen(name);

        %% Compute & Print RMSE
        eV_kf = (permute(xh(3,1:2,:),[2 3 1]) - x(1:2,1:N)).*(permute(xh(3,1:2,:),[2 3 1]) - x(1:2,1:N));
        eV_ekf = (permute(xh(2,1:2,:),[2 3 1]) - x(1:2,1:N)).*(permute(xh(2,1:2,:),[2 3 1]) - x(1:2,1:N));
        eV_ukf = (permute(xh(1,1:2,:),[2 3 1]) - x(1:2,1:N)).*(permute(xh(1,1:2,:),[2 3 1]) - x(1:2,1:N));
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
        pf_err = (permute(xh(4,1:2,:),[2 3 1]) - x(1:2,1:N)).*(permute(xh(4,1:2,:),[2 3 1]) - x(1:2,1:N));
        err = (zV_ukf(:,:) - x(1:2,1:N)).*(zV_ukf(:,:) - x(1:2,1:N));
        RMSE_pf = sqrt(sum(pf_err,2)/T)
        RMSE_meas = sqrt(sum(err,2)/T)
    %     %% Plot results
    %     figure('units','centimeters','position',[.1 .1 10.2 9.1])
    %     %plot(1:N, sqrt(err(:,:)), 'k--', 'LineWidth', 1) 
    %     plot(1:N, pV_err(8,:), 'k--', 'LineWidth', 0.5)
    %     hold on
    %     plot(1:N, pV_err(1,:), 'c-', 'LineWidth', 1)
    %     plot(1:N, pV_err(2,:), 'bx-', 'LineWidth', 1)
    %     plot(1:N, pV_err(3,:), 'r-', 'LineWidth', 1)
    %     plot(1:N, pV_err(4,:), 'g-', 'LineWidth', 1)
    %     plot(1:N, pV_err(5,:), 'mo-', 'LineWidth', 1)
    %     plot(1:N, pV_err(6,:), 'Marker','+','LineStyle','--','Color',[0.8,0.8,0], 'LineWidth', 1)
    %     plot(1:N, pV_err(7,:), 'Marker','x','LineStyle','--','Color',[1,0.6,0.4], 'LineWidth', 1)
    %     h_legend = legend('Meas','KF', 'EKF', 'UKF','PF-CCH','PF-CCHR', 'PF-PCH', 'PF-PCHR', 'Orientation','horizontal');
    %     set(h_legend,'FontSize',9);
    % 
    %     figure('units','centimeters','position',[.1 .1 10.2 9.1])
    %     %plot(1:N, sqrt(err(:,:)), 'k--', 'LineWidth', 1) 
    %     hold on
    %     plot(1:N, hV_err(1,:), 'c-', 'LineWidth', 1)
    %     plot(1:N, hV_err(2,:), 'bx-', 'LineWidth', 1)
    %     plot(1:N, hV_err(3,:), 'r-', 'LineWidth', 1)
    %     plot(1:N, hV_err(4,:), 'g-', 'LineWidth', 1)
    %     plot(1:N, hV_err(5,:), 'mo-', 'LineWidth', 1)
    %     plot(1:N, hV_err(6,:), 'Marker','+','LineStyle','--','Color',[0.8,0.8,0], 'LineWidth', 1)
    %     h_legend = legend('Meas','KF', 'EKF', 'UKF', 'PF-CV','PF-CH','PF-CHR', 'Orientation','horizontal');
    %     h_legend = legend('Meas','KF-CV', 'EKF-CH', 'UKF', 'PF-CV','PF-CH','PF-CHR', 'Orientation','horizontal');
    %     set(h_legend,'FontSize',9);

        PV_err(Run_iter,:,:) = pV_err;
        % for k=1:2                                 % plot results
        % %     figure
        %     subplot(2,1,k)
        %     %figure
        %     if (k<3)
        %         axis([0,N,0,0.30])
        %     elseif (k==3)
        %         axis([0,N,0,0.06])
        %     else
        %         axis([0,N,0,0.1])
        %     end
        %      hold on
        %     plot(1:N, sqrt(err(k,:)), 'k--', 'LineWidth', 1) 
        %     plot(1:N, sqrt(eV_kf(k,:)), 'c-', 'LineWidth', 1)
        %     plot(1:N, sqrt(eV_ekf(k,:)), 'b-', 'LineWidth', 1)
        %     plot(1:N, sqrt(eV_ukf(k,:)), 'r-', 'LineWidth', 1)
        %     plot(1:N, sqrt(pf_err(k,:)), 'g-', 'LineWidth', 1) 
        %     %xlabel('label_{subscript}')
        %     str = sprintf(['RMSE vs. Time for x_{' num2str(k) ',k}']);
        %     title(str)
        %     h_legend = legend('Meas', 'KF', 'EKF', 'UKF', 'PF','Orientation','horizontal');
        %     set(h_legend,'FontSize',9);
        %     %text(0,0.90,'$\textcircled{a}$', 'Interpreter', 'latex');
        % 
        % 
        % end
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

    end
    err = sqrt((sum(PV_err.^2,1))/No_of_Runs);

    %% Plot results
    % Plot error 
    figure('units','centimeters','position',[.1 .1 30 20])
    %plot(1:N, sqrt(err(:,:)), 'k--', 'LineWidth', 1) 
    %plot(1:N, permute(err(1,end,:),[2 3 1]), 'k--', 'LineWidth', 0.5)
    hold on
    h1=plot(1:N, permute(err(1,1,:),[2 3 1]), 'c-', 'LineWidth', 1);
    h2=plot(1:N, permute(err(1,6,:),[2 3 1]), 'g-', 'LineWidth', 1);
    h3=plot(1:N, permute(err(1,2,:),[2 3 1]), 'bx-', 'LineWidth', 1);
    h4=plot(1:N, permute(err(1,3,:),[2 3 1]), 'r-', 'LineWidth', 1);
    %plot(1:N, permute(err(1,4,:),[2 3 1]), 'Marker','+','LineStyle','--','Color',[1,0.6,0], 'LineWidth', 1)
    h5=plot(1:N, permute(err(1,5,:),[2 3 1]), 'mo-', 'LineWidth', 1);
    %plot(1:N, permute(err(1,7,:),[2 3 1]), 'Marker','+','LineStyle','--','Color',[0.8,0.8,0], 'LineWidth', 1)
    % plot(1:N, permute(err(1,7,:),[2 3 1]), 'Marker','x','LineStyle','--','Color',[1,0.6,0.4], 'LineWidth', 1)
    %h_legend = legend('Meas','KF', 'EKF-CCH', 'UKF-CCH', 'PF-CCV', 'PF-CCH','PF-CCV-OPT', 'PF-CCH-OPT', 'Orientation','horizontal');%'PF-PCH', 'PF-PCHR');
    plot(repmat(8, 2), [-5, 5], 'k-');
    plot(repmat(73, 2), [-5, 5], 'k-');
    plot(repmat(93, 2), [-5, 5], 'k-');
    h_legend = legend([h1, h2, h3, h4, h5], {'KF-CV', 'PF-CV', 'EKF-CH', 'UKF-CH', 'PF-CH'}, 'Orientation','horizontal');
    set(h_legend,'FontSize',9);
    tit = "RMS Positional error (\sigma_{m}=";
    tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
    title(tit)
    ylabel("RMSE (m)", 'FontSize', 20, 'FontWeight', 'bold')
    xlabel("Time index t", 'FontSize', 20, 'FontWeight', 'bold')
    set(gca,'FontSize',20)
    set(h_legend,'FontSize',16);
    axis([0 130 0 3])
    box on

    figure('units','centimeters','position',[.1 .1 30 20])
    imagesc([min_x max_x], [min_y max_y], flipud(img));
    hold on
    rectangle('Position',[8.75 4 0.75 3], 'FaceColor', [0.8 0.8 0.8])
    text(8.8,6.8,'A_3', 'FontWeight', 'bold')
    rectangle('Position',[4 6 2 1], 'FaceColor', [0.8 0.8 0.8])
    text(4.05,6.8,'A_2', 'FontWeight', 'bold')
    rectangle('Position',[1 3 2 2], 'FaceColor', [0.8 0.8 0.8])
    text(1.05,4.8,'A_1', 'FontWeight', 'bold')
    h1 = plot(x1(:,1), y1(:,1), 'rx--', 'LineWidth', 2)
    h2 = plot(x_obstructed(:,1), y_obstructed(:,1), 'r--', 'LineWidth', 2)
    h3 = plot(x1(:,2), y1(:,2), 'bx--', 'LineWidth', 2)
    h4 = plot(x_obstructed(:,2), y_obstructed(:,2), 'b--', 'LineWidth', 2)
    h5 = plot(x1(:,3), y1(:,3), 'cx--', 'LineWidth', 2)
    h6 = plot(x_obstructed(:,3), y_obstructed(:,3), 'c--', 'LineWidth', 2)
    set(gca, 'ydir','normal');
    h_legend = legend([h1, h3, h5], {'Target 1', 'Target 2', 'Target3'}, 'Orientation','horizontal');
    %tit = "RMS Positional error (\sigma_{m}=";
    %tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
    %title(tit)
    ylabel("y-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
    xlabel("x-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
    set(gca,'FontSize',20)
    set(h_legend,'FontSize',16);
    axis([0 10 0 10])
    box on

    % Plot trajectory
%     figure('units','centimeters','position',[.1 .1 30 20])
%     hold on
%     h1 = plot(sV_ukf(1,1:k),sV_ukf(2,1:k),'k--','LineWidth',4);
%     h2 = plot(sV_ukf(1,k),sV_ukf(2,k),'ko','MarkerSize', 20, 'LineWidth', 4, 'MarkerFaceColor','red');
%     h2 = plot(sV_ukf(1,1),sV_ukf(2,1),'ko','MarkerSize', 20, 'LineWidth', 4, 'MarkerFaceColor','green');h3 = plot(zV_ukf(1,1:k), zV_ukf(2,1:k),'r.','MarkerSize', 15);
%     h_legend = legend([h1, h3], {'Ground Truth', 'Measurements'}, 'Orientation','Vertical');
%     title("True trajectory of target in 2-D (\sigma_m=1m)")
%     ylabel("y-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
%     xlabel("x-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
%     set(gca,'FontSize',20)
%     set(h_legend,'FontSize',16);
%     box on
    
    RMSE_mean(:, r_iter) = sum(permute(err(1,:,:),[2 3 1]),2)/N;
    RMSE_perc(:, r_iter) = RMSE_mean(:, r_iter)/r_list(r_iter);
end