r_list = [0.5, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
Sim_No = 10;
clear F
clear pf
PV_err = [];
PV_lost = [];
RMSE_test = [];
lost_test = [];
for test_iter = 1:size(r_list,2)
    RMSE_Sim = [];
    lost_Sim = [];
    for sim_iter = 1:Sim_No
        close all
        RMSE_avg = [];
        lost_avg = [];
        % Number of Target Tracks
        TrackNum = 1 ;
        [DataList,x1,y1] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, r_list(1,test_iter), 16, 0);

        Dt=1;
%         q = 0.01;
        % Initiate KF parameters
        nx=4;      %number of state dims
        ny = 2;    %number of observation dims
        n=4;      %number of state
        %q=0.03;    %std of process 
        r=r_list(1,test_iter);    %std of measurement
        %s.Q=[1^3/3, 0, 1^2/2, 0;  0, 1^3/3, 0, 1^2/2; 1^2/2, 0, 1, 0; 0, 1^2/2, 0, 1]*q^2; % covariance of process
        s.R=r^2*eye(n/2);        % covariance of measurement  
%         s.sys=(@(x)[x(1)+ x(3); x(2)+x(4); x(3); x(4)]);  % assuming measurements arrive 1 per sec
%         s.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2; % covariance of process
        s.sys=@(x)[x(1)+ x(3)*cos(x(4)); x(2)+x(3)*sin(x(4)); x(3); x(4)];  % nonlinear state equations
        s.Q = diag([0,0,0.01^2,0.3^2]);
        s.obs=@(x)[x(1);x(2)];                               % measurement equation
        st=[x_true(1,1);y_true(1,1)];                                % initial state
        s.x_init = [];
        s.P_init = [];
        for i = 1:TrackNum
            % s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); 0; DataList(2,i,2)-DataList(2,i,1)]; %initial state
            s.x_init(:,i)=[x1(2,1); y1(2,1); sqrt((x1(2,1)-x1(1,1))^2+(y1(2,1)-y1(1,1))^2); atan((y1(2,1)+0.1-y1(1,1))/(x1(2,1)+0.1-x1(1,1)))]; %initial state
        end
        s.P_init = diag([r^2, r^2, 2*r^2, 2*r^2]);                               % initial state covraiance

        % Process and Observation handles and covariances
        TrackObj.sys          = s.sys;
        TrackObj.obs          = s.obs;
        TrackObj.Q          = s.Q;
        TrackObj.R          = s.R;

        % initial covariance assumed same for all tracks
        TrackObj.P          = s.P_init;

        N=size(DataList,2) ;                                    % total dynamic steps

        % Instantiate EKF and vars to store output
        xV_ekf = zeros(n,N);          %estmate        % allocate memory
        xV_ekf(:,1) = s.x_init(:,1);
        PV_ekf = zeros(1,N);    
        %PV_ekf = cell(1,N);     % use to display ellipses 
        sV_ekf = zeros(n/2,N);          %actual
        sV_ekf(:,1) = st;
        zV_ekf = zeros(n/2,N);
        eV_ekf = zeros(n/2,N);

         %% Initiate PF parameters

        % Process equation x[k] = sys(k, x[k-1], u[k]);
        sys_cch = @(k, xkm1, uk) [xkm1(1)+1*xkm1(3)*cos(xkm1(4)); xkm1(2)+1*xkm1(3)*sin(xkm1(4)); xkm1(3)+ uk(3); xkm1(4) + uk(4)];

        % Observation equation y[k] = obs(k, x[k], v[k]);
        obs = @(k, xk, vk) [xk(1)+vk(1); xk(2)+vk(2)];                  % (returns column vector)

        % PDF of process and observation noise generator function
        nu = 4;                                           % size of the vector of process noise
        %sigma_u = q;
        %cov_u = [Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, 1]*sigma_u^2;
        gen_sys_noise_cch = @(u) mvnrnd(zeros(1, nu), diag([0,0,0.01^2,0.3^2])); 
        % PDF of observation noise and noise generator function
        nv = 2;                                           % size of the vector of observation noise
        sigma_v = r;
        cov_v = sigma_v^2*eye(nv);
        p_obs_noise   = @(v) mvnpdf(v, zeros(1, nv), cov_v);
        gen_obs_noise = @(v) mvnrnd(zeros(1, nv), cov_v);         % sample from p_obs_noise (returns column vector)

        % Initial PDF
        gen_x0_cch = @(x) mvnrnd([obs_x(1,1), obs_y(1,1), 0, 0],diag([100^2, 100^2, 0, 0]));

        % Observation likelihood PDF p(y[k] | x[k])
        % (under the suposition of additive process noise)
        p_yk_given_xk = @(k, yk, xk) p_obs_noise((yk - obs(k, xk, zeros(1, nv)))');

        % Separate memory space
        %x = [x_true(:,1), y_true(:,1)]';
        %y = [obs_x(:,1)'; obs_y(:,1)']; % True state and observations

        % Assign PF parameter values
        pf.k               = 1;                   % initial iteration number
        pf.Np              = 1000;                 % number of particles
        %pf.w               = zeros(pf.Np, T);     % weights
        pf.particles       = zeros(5, pf.Np); % particles
        %pf.p_x0 = p_x0;                          % initial prior PDF p(x[0])
        %pf.p_xk_given_ xkm1 = p_xk_given_xkm1;   % transition prior PDF p(x[k] | x[k-1])
        %pf.xhk = s.x;
        pf.resampling_strategy = 'systematic_resampling';

        % PF-PCHR
        %pf_pchr = ParticleFilterMin(pf);

        % PF-CCH
        pf.sys = sys_cch;
        pf.particles = zeros(nx, pf.Np); % particles
        pf.gen_x0 = gen_x0_cch;
        pf.obs = p_yk_given_xk;
        pf.obs_model = @(xk) [xk(1,:); xk(2,:)];
        pf.R = cov_v;
        pf.clutter_flag = 1;
        pf.sys_noise = gen_sys_noise_cch;
        %s.Q = diag([0.01^10,0.01^10,0.01^2,0.3^2]);
        %pf.kf = UKalmanFilter(s, 0.5, 0, 2);

        % Initiate Tracklist
        TrackList = [];
        TrackListPF = [];

        %% Estimated State container PF
        xh = zeros(TrackNum,nx, N);
        RMSE = zeros(TrackNum,N);
        lost = zeros(TrackNum,N);
        for i=1:TrackNum,

            TrackObj.x          = s.x_init(:,i);
            TrackList{i}.TrackObj = TrackObj;
            pf.gen_x0 = @(x) mvnrnd([x1(2,1); y1(2,1); sqrt((x1(2,1)-x1(1,1))^2+(y1(2,1)-y1(1,1))^2); atan((y1(2,1)+0.1-y1(1,1))/(x1(2,1)+0.1-x1(1,1)))],diag([sigma_v^2, sigma_v^2, 2*sigma_v^2, 2*sigma_v^2]));
            pf.xhk = [s.x_init(1,i),s.x_init(2,i),0,0]';
            TrackListPF{i}.TrackObj = ParticleFilterMin(pf);
            xh(i,:,1) = pf.xhk;
            %TrackList{i}.TrackObj.x(ObservInd) = CenterData(:,i);

        end;

        ekf = EKalmanFilter(TrackObj);

        img = imread('maze.png');

        % set the range of the axes
        % The image will be stretched to this.
        min_x = 0;
        max_x = 10;
        min_y = 0;
        max_y = 10;

        % make data to plot - just a line.
        %x = min_x:max_x;
        %y = (6/8)*x;

        figure('pos',[10 10 1610 810])
        for i=2:N
            % Compute squared error
            %eV_ekf(:,i) = (TrackList{1,1}.TrackObj.x(1:2,1) - st).*(TrackList{1,1}.TrackObj.x(1:2,1) - st);
            clf; % Flip the image upside down before showing it
            %imagesc([min_x max_x], [min_y max_y], flipud(img));
            axis([0 25 0 20])
            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            
  
            
            fprintf('R = %d\n',r_list(1,test_iter));
            fprintf('\nSim_No = %d/%d\n',sim_iter,Sim_No);
            fprintf('Iteration = %d/%d\n',i,N);
            [TrackList, Validation_matrix] = Observation_Association(TrackList, DataList{i}, ekf);
            TrackList = PDAF_Update(TrackList, DataList{i}, Validation_matrix);

            for t = 1:TrackNum
                TrackListPF{t}.TrackObj.pf.k = i;
                TrackListPF{t}.TrackObj.pf.z = DataList{i};
                TrackListPF{t}.TrackObj.pf = TrackListPF{t}.TrackObj.Predict(TrackListPF{t}.TrackObj.pf);
                TrackListPF{t}.TrackObj.pf = TrackListPF{t}.TrackObj.Update(TrackListPF{t}.TrackObj.pf);
                xh(t,:,i) = TrackListPF{t}.TrackObj.pf.xhk;
                RMSE(t,i) = sqrt((xh(t,1,i) - x1(i,1))^2 + (xh(t,2,i)-y1(i,1))^2);
            end

            xV_ekf(:,i) = TrackList{1,1}.TrackObj.x;
            RMSE(TrackNum+1,i) = sqrt((xV_ekf(1,i) - x1(i,1))^2 + (xV_ekf(2,i)-y1(i,1))^2);
            for l = 1:size(RMSE,1)
                lost(l,i) = RMSE(l,i) > 3;
            end
            st = [x_true(i,1); y_true(i,1)];
            sV_ekf(:,i)= st;
            h1 = plot(x1(1:i,1),y1(1:i,1),'b.-','LineWidth',1);
            h4 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
            h6 = plot(DataList{i}(1,1),DataList{i}(2,1),'r*','MarkerSize', 10);
            h2 = plot(permute(xh(1,1,1:i),[2 3 1]),permute(xh(1,2,1:i),[2 3 1]),'c.-','LineWidth',1);
            h3 = plot(xV_ekf(1,1:i),xV_ekf(2,1:i),'r.-','LineWidth',1);
            h5 = plot(TrackListPF{1}.TrackObj.pf.particles(1,:),TrackListPF{1}.TrackObj.pf.particles(2,:),'k.', 'MarkerSize', 1);
            h10 = plot_gaussian_ellipsoid([permute(xh(1,1,i),[2 3 1]);permute(xh(1,2,i),[2 3 1])], [std(TrackListPF{1}.TrackObj.pf.particles(1,:))^2,0;0,std(TrackListPF{1}.TrackObj.pf.particles(2,:))^2]); 
            h7 = plot(x1(i,1),y1(i,1),'bo','MarkerSize', 20);
            h8 = plot(permute(xh(1,1,i),[2 3 1]),permute(xh(1,2,i),[2 3 1]),'co','MarkerSize', 20, 'LineWidth',2);
            h9 = plot(xV_ekf(1,i),xV_ekf(2,i),'ro','MarkerSize', 20, 'LineWidth',2);
            h10 = plot_gaussian_ellipsoid([xV_ekf(1,i);xV_ekf(2,i)], TrackList{1,1}.TrackObj.P(1:2,1:2));
            h_legend = legend([h1 h4 h3 h2],'Ground Truth', 'measurements', 'PDAF-EKF', 'PDAF-PF');%, 'PF-PCH', 'PF-PCHR');
            set(h_legend, 'FontSize', 18, 'Orientation','horizontal');
            %      % set the y-axis back to normal.
            %tw = sprintf('True vs Estimated Trajectory(?_{m}=%1.1f, ?_{c}=%d',[r,50]);
            title(sprintf('True vs Estimated Trajectory (\\sigma_{m}=%1.1f m, \\lambda_{c}=%d)',r,50),'fontsize',18);
            xlabel('X-Coordinate (m)','fontsize',18);
            ylabel('Y-Coordinate (m)','fontsize',18); 
            %set(gca,'ydir','normal');
            %axis([0 10 0 10])
            pause(0.1)
%             ax = gca;
%             ax.Units = 'pixels';
%             pos = ax.Position;
%             marg = 30;
%             rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
%             F(i) = getframe(gcf);
%             ax.Units = 'normalized';
            
        end
%         F = F(2:end);
%         vidObj = VideoWriter('myPeaks4.avi');
%         vidObj.Quality = 100;
%         vidObj.FrameRate = 10;
%         open(vidObj);
%         writeVideo(vidObj, F);
%         close(vidObj);
%         winopen('myPeaks4.avi');
        RMSE_Sim(:,end+1) = mean(RMSE,2);
        lost_Sim(:,end+1) = mean(lost,2);
        PV_err(sim_iter,:,:) = RMSE;
        PV_lost(sim_iter,:,:) = lost;
%         RMSE_avg(:,end+1) = mean(RMSE,2);
%         lost_avg(:,end+1) = mean(lost,2);
    end
    err = sqrt((sum(PV_err.^2,1))/Sim_No);
    lost_p = sqrt((sum(PV_lost.^2,1))/Sim_No);
    RMSE_test(:,end+1) = mean(RMSE_Sim,2);
    lost_test(:,end+1) = mean(lost_Sim,2);
end
% Compute & Print RMSE
% RMSE_avg(:,end+1) = mean(RMSE,2)

figure
hold on;
h2 = plot(sV_ekf(1,1:N),sV_ekf(2,1:N),'b.-','LineWidth',1);
h4 = plot(xV_ekf(1,1:N),xV_ekf(2,1:N),'r.-','LineWidth',1);
