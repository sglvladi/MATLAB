r_list = [1];
lambda_list = 1;
Sim_No = 10;
clear F
PV_err = [];
PV_lost = [];
RMSE_test = [];
lost_test = [];
for test_iter = 1:size(lambda_list,2)
    RMSE_Sim = [];
    lost_Sim = [];
    for sim_iter = 1:Sim_No
        close all
        RMSE_avg = [];
        lost_avg = [];
        % Number of Target Tracks
        TrackNum = 1 ;
        [DataList,x1,y1] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, r_list(test_iter), 16, lambda_list(1,test_iter),2);

        Dt=1;
        q = 0.01;
        % Initiate KF parameters
        nx=4;      %number of state dims
        ny = 2;    %number of observation dims
        n=4;      %number of state
        %q=0.03;    %std of process 
        r=r_list(test_iter);    %std of measurement
        %s.Q=[1^3/3, 0, 1^2/2, 0;  0, 1^3/3, 0, 1^2/2; 1^2/2, 0, 1, 0; 0, 1^2/2, 0, 1]*q^2; % covariance of process
        s.R=r^2*eye(n/2);        % covariance of measurement  
%         s.sys=(@(x)[x(1)+ x(3); x(2)+x(4); x(3); x(4)]);  % assuming measurements arrive 1 per sec
%         s.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2; % covariance of process
        s.sys=@(x)[x(1)+ x(3)*cos(x(4)); x(2)+x(3)*sin(x(4)); x(3); x(4)];  % nonlinear state equations
        s.Q = diag([0.01,0.01,0.01^2,0.3^2]);
        s.obs=@(x)[x(1);x(2)];                               % measurement equation
        st=[x_true(1,1);y_true(1,1)];                                % initial state
        s.x_init = [];
        s.P_init = [];
        for i = 1:TrackNum
            % s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); 0; DataList(2,i,2)-DataList(2,i,1)]; %initial state
            s.x_init(:,i)=[x1(2,1); y1(2,1); sqrt((x1(2,1)-x1(1,1))^2+(y1(2,1)-y1(1,1))^2); atan((y1(2,1)-y1(1,1))/(x1(2,1)-x1(1,1)))-pi]; %initial state
        end
        s.P_init = diag([q^2, q^2, 2*q^2, 0.1]);                               % initial state covraiance

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

        %% Initiate PF parameters
        nx = 4;      % number of state dims
        nu = 4;      % size of the vector of process noise
        nv = 2;      % size of the vector of observation noise
        % Prior PDF generator
        gen_x0_cch = @(Np) mvnrnd(repmat([0,0,0,0],Np,1),diag([q^2, q^2, 2*q^2, 0.1]));
        % Process equation x[k] = sys(k, x[k-1], u[k]);
        sys_cch = @(k, xkm1, uk) [xkm1(1,:)+1*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+1*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)'];
        % PDF of process noise generator function
        gen_sys_noise_cch = @(u) mvnrnd(zeros(size(u,2), nu), diag([0,0,q^2,0.3^2])); 
        % Observation equation y[k] = obs(k, x[k], v[k]);
        obs = @(k, xk, vk) [xk(1)+vk(1); xk(2)+vk(2)];                  % (returns column vector)
        % PDF of observation noise and noise generator function
        sigma_u = q;
        sigma_v = r;
        cov_v = sigma_v^2*eye(nv);
        p_obs_noise   = @(v) mvnpdf(v, zeros(1, nv), cov_v);
        gen_obs_noise = @(v) mvnrnd(zeros(1, nv), cov_v);         % sample from p_obs_noise (returns column vector)
        % Observation likelihood PDF p(y[k] | x[k])
        % (under the suposition of additive process noise)
        p_yk_given_xk = @(k, yk, xk) p_obs_noise((yk - obs(k, xk, zeros(1, nv)))');
        % Assign PF parameter values
        pf.k               = 1;                   % initial iteration number
        pf.Np              = 10000;                 % number of particles
        pf.particles       = zeros(5, pf.Np); % particles
        pf.resampling_strategy = 'systematic_resampling';
        pf.sys = sys_cch;
        pf.particles = zeros(nx, pf.Np); % particles
        pf.gen_x0 = gen_x0_cch(pf.Np);
        pf.obs = p_yk_given_xk;
        pf.obs_model = @(xk) [xk(1,:); xk(2,:)];
        pf.R = cov_v;
        pf.clutter_flag = 1;
        %pf.multi_flag = 1;
        pf.sys_noise = gen_sys_noise_cch;

        % PF-PCHR
        %pf_pchr = ParticleFilterMin(pf);

        % PF-CCH
%         pf.sys = sys_cch;
%         pf.particles = zeros(nx, pf.Np); % particles
%         pf.gen_x0 = gen_x0_cch;
%         pf.obs = p_yk_given_xk;
%         pf.obs_model = @(xk) [xk(1,:); xk(2,:)];
%         pf.R = cov_v;
%         pf.clutter_flag = 1;
%         pf.sys_noise = gen_sys_noise_cch;

        % Initiate Tracklist
        TrackList = [];
        TrackListKF = [];
        TrackListUKF = [];
        TrackListPF = [];

        %% Estimated State container PF
        xh = zeros(TrackNum,nx, N);
        RMSE = zeros(TrackNum,N);
        lost = zeros(TrackNum,N);
        
        ekf = EKalmanFilter(TrackObj);
        TrackObj.sys=@(x,u)[x(1)+ x(3)*cos(x(4))+u(1); x(2)+x(3)*sin(x(4))+u(2); x(3)+u(3); x(4)+u(4)];
        TrackObj.obs=@(x,u)[x(1)+u(1);x(2)+u(2)];
        ukf = UKalmanFilter(TrackObj, 0.5, 0, 2);
        TrackObj.Q=[Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt]*q^2; % covariance of process
        TrackObj.sys=[1 0 Dt 0; 0 1 0 Dt; 0 0 1 0; 0 0 0 1];
        TrackObj.obs = [1 0 0 0; 0 1 0 0];
        TrackObj.x(:,i)= [x_true(2); y_true(2); (x_true(2)-x_true(1)); (y_true(2)-y_true(1))]; %initial state
        TrackObj.P = diag([q^2, q^2, 2*q^2, 2*q^2]);
        kf = KalmanFilter_new(TrackObj);
        for i=1:TrackNum

            TrackObj.x          = s.x_init(:,i);
            TrackList{i}.TrackObj = ekf.s;
            TrackListKF{i}.TrackObj = kf.s;
            TrackListUKF{i}.TrackObj = ukf.s;
            pf.gen_x0 = @(x) mvnrnd(repmat([x1(2,1); y1(2,1); sqrt((x1(2,1)-x1(1,1))^2+(y1(2,1)-y1(1,1))^2); atan((y1(2,1)-y1(1,1))/(x1(2,1)-x1(1,1)))-pi]',pf.Np,1),diag([sigma_u^2, sigma_u^2, 2*sigma_u^2, 0.1]));
            pf.xhk = [s.x_init(1,i),s.x_init(2,i),0,0]';
            TrackListPF{i}.TrackObj = ParticleFilterMin2(pf);
            xh(i,:,1) = pf.xhk;
            %TrackList{i}.TrackObj.x(ObservInd) = CenterData(:,i);

        end;

        
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
            clf; % Flip the image upside down before showing it
            %imagesc([min_x max_x], [min_y max_y], flipud(img));
            %axis([0 15 0 10])
            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            
            fprintf('\nLambda = %d\n',lambda_list(1,test_iter));
            fprintf('Sim_No = %d/%d\n',sim_iter,Sim_No);
            fprintf('Iteration = %d/%d\n',i,N);
            
            [TrackListKF, Validation_matrix] = Observation_Association(TrackListKF, DataList{i}, kf);
            TrackListKF = PDAF_Update(TrackListKF, DataList{i}, Validation_matrix);
            
            [TrackList, Validation_matrix] = Observation_Association(TrackList, DataList{i}, ekf);
            TrackList = PDAF_Update(TrackList, DataList{i}, Validation_matrix);
            
            [TrackListUKF, Validation_matrix] = Observation_Association(TrackListUKF, DataList{i}, ukf);
            TrackListUKF = PDAF_UKF_Update(TrackListUKF, DataList{i}, Validation_matrix);

            for t = 1:TrackNum
                TrackListPF{t}.TrackObj.pf.k = i;
                TrackListPF{t}.TrackObj.pf.z = DataList{i};
                TrackListPF{t}.TrackObj.pf = TrackListPF{t}.TrackObj.Predict(TrackListPF{t}.TrackObj.pf);
                TrackListPF{t}.TrackObj.pf = TrackListPF{t}.TrackObj.Update(TrackListPF{t}.TrackObj.pf);
                xh(t,:,i) = TrackListPF{t}.TrackObj.pf.xhk;
                RMSE(t,i) = sqrt((xh(t,1,i) - x1(i,1))^2 + (xh(t,2,i)-y1(i,1))^2);
            end

            xV_kf(:,i) = TrackListKF{1,1}.TrackObj.x;
            xV_ekf(:,i) = TrackList{1,1}.TrackObj.x;
            xV_ukf(:,i) = TrackListUKF{1,1}.TrackObj.x;
            RMSE(TrackNum+1,i) = sqrt((xV_ekf(1,i) - x1(i,1))^2 + (xV_ekf(2,i)-y1(i,1))^2);
            RMSE(TrackNum+2,i) = sqrt((xV_ukf(1,i) - x1(i,1))^2 + (xV_ukf(2,i)-y1(i,1))^2);
            RMSE(TrackNum+3,i) = sqrt((xV_kf(1,i) - x1(i,1))^2 + (xV_kf(2,i)-y1(i,1))^2);
%             for l = 1:size(RMSE,1)
%                 lost(l,i) = RMSE(l,i) > 3*r_list(test_iter);
%             end
            lost = RMSE > 3*r_list(test_iter);
            st = [x_true(i,1); y_true(i,1)];
            sV_ekf(:,i)= st;
            % Compute squared error
            %eV_ekf(:,i) = (TrackList{1,1}.TrackObj.x(1:2,1) - st).*(TrackList{1,1}.TrackObj.x(1:2,1) - st);
            
%             h1 = plot(x1(1:i,1),y1(1:i,1),'k.-','LineWidth',1);
%             h2 = plot(permute(xh(1,1,1:i),[2 3 1]),permute(xh(1,2,1:i),[2 3 1]),'m.-','LineWidth',1);
%             h3 = plot(xV_kf(1,1:i),xV_kf(2,1:i),'c.-','LineWidth',1);
%             h4 = plot(xV_ekf(1,1:i),xV_ekf(2,1:i),'b.-','LineWidth',1);
%             h5 = plot(xV_ukf(1,1:i),xV_ukf(2,1:i),'r.-','LineWidth',1);
%             h6 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
%             h7 = plot(TrackListPF{1}.TrackObj.pf.particles(1,:),TrackListPF{1}.TrackObj.pf.particles(2,:),'k.', 'MarkerSize', 1);
%             h8 = plot(DataList{i}(1,1),DataList{i}(2,1),'r*','MarkerSize', 10);
%             h9 = plot(x1(i,1),y1(i,1),'ko','MarkerSize', 20);
%             h10 = plot(permute(xh(1,1,i),[2 3 1]),permute(xh(1,2,i),[2 3 1]),'mo','MarkerSize', 20, 'LineWidth',2);
%             h11 = plot(xV_kf(1,i),xV_kf(2,i),'co','MarkerSize', 20, 'LineWidth',2);
%             h12 = plot(xV_ekf(1,i),xV_ekf(2,i),'bo','MarkerSize', 20, 'LineWidth',2);
%             h13 = plot(xV_ukf(1,i),xV_ukf(2,i),'ro','MarkerSize', 20, 'LineWidth',2);
%             h14 = plot_gaussian_ellipsoid([xV_ekf(1,i);xV_ekf(2,i)], TrackList{1,1}.TrackObj.P(1:2,1:2));
%             h_legend = legend([h1 h5 h3 h4 h5 h2],'Ground Truth', 'measurements', 'PDAF-KF', 'PDAF-EKF', 'PDAF-UKF', 'PDAF-PF');%, 'PF-PCH', 'PF-PCHR');
%             set(h_legend, 'FontSize', 18, 'Orientation','horizontal');
%             %      % set the y-axis back to normal.
%             %tw = sprintf('True vs Estimated Trajectory(?_{m}=%1.1f, ?_{c}=%d',[r,50]);
%             title(sprintf('True vs Estimated Trajectory (\\sigma_{m}=%1.1f m, \\lambda_{c}=%d)',r,50),'fontsize',18);
%             xlabel('X-Coordinate (m)','fontsize',18);
%             ylabel('Y-Coordinate (m)','fontsize',18); 
%             %set(gca,'ydir','normal');
%             axis([-5 30 -5 30])
%              pause(0.001)
%             ax = gca;
%             ax.Units = 'pixels';
%             pos = ax.Position;
%             marg = 30;
%             rect = [-marg, -marg, pos(3)+2*marg, pos(4)+2*marg];
%             F(i) = getframe(gcf);
%             ax.Units = 'normalized';
            
        end
%         F = F(2:end);
%         vidObj = VideoWriter(sprintf('test%d.avi',sim_iter));
%         vidObj.Quality = 100;
%         vidObj.FrameRate = 10;
%         open(vidObj);
%         writeVideo(vidObj, F);
%         close(vidObj);
%         winopen(sprintf('test%d.avi',sim_iter));
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

% Plot error 
    %N=124
    %err = err(:,:, 1:3:end)
    figure('units','centimeters','position',[.1 .1 30 20])
    %plot(1:N, sqrt(err(:,:)), 'k--', 'LineWidth', 1) 
    %plot(1:N, permute(err(1,end,:),[2 3 1]), 'k--', 'LineWidth', 0.5)
    hold on
    h1=plot(1:N, permute(RMSE(1,:),[2 3 1]), 'mo-', 'LineWidth', 1);
    h2=plot(1:N, permute(RMSE(3,:),[2 3 1]), 'r-', 'LineWidth', 1);
    h3=plot(1:N, permute(RMSE(2,:),[2 3 1]), 'bx-', 'LineWidth', 1);
    %h4=plot(1:N, permute(err(1,3,:),[2 3 1]), 'r-', 'LineWidth', 1);
    %plot(1:N, permute(err(1,4,:),[2 3 1]), 'Marker','+','LineStyle','--','Color',[1,0.6,0], 'LineWidth', 1)
    %h5=plot(1:N, permute(err(1,5,:),[2 3 1]), 'mo-', 'LineWidth', 1);
    %plot(1:N, permute(err(1,7,:),[2 3 1]), 'Marker','+','LineStyle','--','Color',[0.8,0.8,0], 'LineWidth', 1)
    % plot(1:N, permute(err(1,7,:),[2 3 1]), 'Marker','x','LineStyle','--','Color',[1,0.6,0.4], 'LineWidth', 1)
    %h_legend = legend('Meas','KF', 'EKF-CCH', 'UKF-CCH', 'PF-CCV', 'PF-CCH','PF-CCV-OPT', 'PF-CCH-OPT', 'Orientation','horizontal');%'PF-PCH', 'PF-PCHR');
    plot(repmat(8, 2), [-5, 5], 'k-');
    plot(repmat(73, 2), [-5, 5], 'k-');
    plot(repmat(93, 2), [-5, 5], 'k-');
    h_legend = legend([h2, h3, h1], {'PDAF-EKF-CH', 'PDAF-UKF-CH', 'MC-PDAF-CV'}, 'Orientation','horizontal');
    set(h_legend,'FontSize',9);
    tit = "RMS Positional error (\sigma_{m}=";
    tit = strcat(tit, sprintf("%0.2fm)",r_list(test_iter)));
    title(tit)
    ylabel("RMSE (m)", 'FontSize', 20, 'FontWeight', 'bold')
    xlabel("Time index t", 'FontSize', 20, 'FontWeight', 'bold')
    set(gca,'FontSize',20)
    set(h_legend,'FontSize',16);
    axis([0 130 0 3])
    box on

    
figure('units','centimeters','position',[.1 .1 30 20])
hold on
h1 = plot(x_true(1:100),y_true(1:100),'k--','LineWidth',4);
h2 = plot(x_true(100),y_true(100),'ko','MarkerSize',20);
h3 = plot(DataList{100}(1,:),DataList{100}(2,:),'r*','MarkerSize', 10);
h4 = plot(DataList{100}(1,1),DataList{100}(2,1),'g*','MarkerSize', 10);
h_legend = legend([h1,h4,h3], {'Ground trouth','True Measurement', 'False Alarms'}); 
set(h_legend,'FontSize',16);
tit = "Simulation snapshot at t=100s";
%tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
title(tit)
ylabel("y-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
xlabel("x-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
set(gca,'FontSize',20)
box on

lost_test=sort(lost_test,2);
figure('units','centimeters','position',[.1 .1 30 20])
hold on
h1 = plot([0, r_list], [0, lost_test(1,:)]*100, 'mo-')
h2 = plot([0, r_list], [0, lost_test(2,:)]*100, 'bx-')
h3 = plot([0, r_list], [0, lost_test(3,:)]*100, 'r-')
h4 = plot([0, r_list], [0, lost_test(4,:)]*100, 'c-')
h_legend = legend([h4,h2,h3, h1], {'PDAF-KF-CV','PDAF-EKF-CH', 'PDAF-UKF-CH', 'MCPDAF-CH'}, 'Orientation', 'Horizontal'); 
set(h_legend,'FontSize',16);
tit = "Probability of track loss vs Measurement noise (\lambda_{FA}V=50)";
%tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
title(tit)
ylabel("Track loss probability (%)", 'FontSize', 20, 'FontWeight', 'bold')
xlabel("\sigma_m (m)", 'FontSize', 20, 'FontWeight', 'bold')
set(gca,'FontSize',20)
axis([0 1 0 100]);
box on

figure('units','centimeters','position',[.1 .1 30 20])
hold on
h1 = plot(lambda_list, lost_test(1,:)*100, 'mo-')
h2 = plot(lambda_list, lost_test(2,:)*100, 'bx-')
h3 = plot(lambda_list, lost_test(3,:)*100, 'r-')
h4 = plot(lambda_list, lost_test(4,:)*100, 'c-')
h_legend = legend([h4,h2,h3, h1], {'PDAF-KF-CV','PDAF-EKF-CH', 'PDAF-UKF-CH', 'MCPDAF-CH'}, 'Orientation', 'Horizontal'); 
set(h_legend,'FontSize',16);
tit = "Probability of track loss vs Clutter rate (\sigma_m=0.5)";
%tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
title(tit)
ylabel("Track loss probability (%)", 'FontSize', 20, 'FontWeight', 'bold')
xlabel("\lambda_{FA}", 'FontSize', 20, 'FontWeight', 'bold')
set(gca,'FontSize',20)
axis([0 100 0 100]);
box on

figure
hold on;
h2 = plot(sV_ekf(1,1:N),sV_ekf(2,1:N),'b.-','LineWidth',1);
h4 = plot(xV_ekf(1,1:N),xV_ekf(2,1:N),'r.-','LineWidth',1);
