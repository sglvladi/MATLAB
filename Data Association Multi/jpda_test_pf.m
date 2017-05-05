% Number of Target Tracks
TrackNum = 3;

% Generate observations
%[DataList,x1,y1] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, 0.1, 2, 100, 1);
RMSE_ekf = zeros(2, TrackNum);

% Number of simulations
SimNum = 1

% Show plots of data
ShowPlots = 1;
SkipFrames = 0; % Skip n frames between consecutive plots (speed-up)
for sim = 1:SimNum

    % Initiate KF parameters
     n=4;      %number of state
     q=0.01;    %std of process 
     r=0.1;    %std of measurement
     s.Q=[1^3/3, 0, 1^2/2, 0;  0, 1^3/3, 0, 1^2/2; 1^2/2, 0, 1, 0; 0, 1^2/2, 0, 1]*q^2; % covariance of process
     s.R=r^2*eye(n/2);        % covariance of measurement  
     s.sys=(@(x)[x(1)+ x(3); x(2)+x(4); x(3); x(4)]);  % assuming measurements arrive 1 per sec
     s.obs=@(x)[x(1);x(2)];                               % measurement equation
     st=[x1(1,:);y1(1,:)];                                % initial state
%     s.x_init = [0,0,0,0];
     s.P_init = [];
% 
%     % Two-point difference initiation
%     %for i = 1:TrackNum
%     %    s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); DataList(1,i,2)-DataList(1,i,1); DataList(2,i,2)-DataList(2,i,1)]; %initial state
%     %end
%     %s.P_init = [q^2, 0, q^2, 0;
%     %            0, q^2, 0, q^2;
%     %            q^2, 0, 2*q^2, 0;
%     %            0, q^2, 0, 2*q^2];                               % initial state covraiance
% 
%     % Single-point initiation
     Vmax = 0.4; % Max velocity = 0.4 m/s
    for i = 1:TrackNum
        s.x_init(:,i)=[x1(1,i);y1(1,i); 0; 0]; %initial state
    end
     s.P_init = diag([q^2, q^2, (Vmax^2/3), (Vmax^2/3)]);
% 
% 
     % Process and Observation handles and covariances
     TrackObj.sys          = s.sys;
     TrackObj.obs          = s.obs;
     TrackObj.Q          = s.Q;
     TrackObj.R          = s.R;
% 
%     % initial covariMeasInd, TrackIndance assumed same for all tracks
     TrackObj.P          = s.P_init;
% 
     N=size(DataList,2) ;                                    % total dynamic steps
% 
     % Instantiate EKF and vars to store output
     Logs = [];
     Log.xV_ekf = zeros(n,N);          %estmate        % allocate memory
     Log.PV_ekf = zeros(1,N);
     %PV_ekf = cell(1,N);     % use to display ellipses 
     Log.sV_ekf = zeros(n/2,N);          %actual
     Log.zV_ekf = zeros(n/2,N);
     Log.eV_ekf = zeros(n/2,N);
     
    %% Initiate PF parameters
    nx=4;      %number of state dims
    % Process equation x[k] = sys(k, x[k-1], u[k]);
    sys_cch = @(k, xkm1, uk) [xkm1(1,:)+1*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+1*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)'];

    % Observation equation y[k] = obs(k, x[k], v[k]);
    obs = @(k, xk, vk) [xk(1)+vk(1); xk(2)+vk(2)];                  % (returns column vector)

    % PDF of process and observation noise generator function
    nu = 4;                                           % size of the vector of process noise
    %sigma_u = q;
    %cov_u = [Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, 1]*sigma_u^2;
    gen_sys_noise_cch = @(u) mvnrnd(zeros(size(u,2), nu), diag([0,0,q^2,0.3^2])); 
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
    pf.Np              = 5000;                 % number of particles
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
    pf.multi_flag = 1;
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
        pf.gen_x0 = @(Np) mvnrnd(repmat([s.x_init(1,i),s.x_init(2,i),0,0],Np,1),diag([q^2, q^2, (Vmax^2/3), 2*sigma_v^2]));
        pf.xhk = [s.x_init(1,i),s.x_init(2,i),0,0]';
        pf.ExistProb = 0.8;
        TrackListPF{i}.TrackObj = ParticleFilterMin2(pf);
        xh(i,:,1) = pf.xhk;
        %TrackList{i}.TrackObj.x(ObservInd) = CenterData(:,i);

    end; 

    %ekf = EKalmanFilter(TrackObj);
    ekf = UKalmanFilter(TrackObj, 0.5, 1, 2);

    img = imread('maze.png');

    % set the range of the axes
    % The image will be stretched to this.
    min_x = 0;
    max_x = 10;
    min_y = 0;
    max_y = 10;

    % make data to plot - just a line.
    x = min_x:max_x;
    y = (6/8)*x;

    figure
    for i=1:N
        fprintf('Iteration = %d/\n',i);
        tempDataList = DataList{i}(:,:);
        tempDataList( :, ~any(tempDataList,1) ) = []; 
%         [TrackList, ValidationMatrix, bettaNTFA] = Observation_Association(TrackList, tempDataList, ekf);
%         %TrackList = JPDAF_EKF_Update(TrackList, DataList(:,:,i), ValidationMatrix', bettaNTFA);
%         [TrackList, betta] = JPDAF_UKF_Update(TrackList, tempDataList, ValidationMatrix', bettaNTFA);
        %TrackList = Track_InitConfDel(TrackList,tempDataList,ValidationMatrix',bettaNTFA, betta);
        ValidationMatrix = zeros(TrackNum, size(tempDataList,2)); 
        tot_gate_area = 0;
        for t = 1:TrackNum
            TrackListPF{t}.TrackObj.pf.k = i;
            TrackListPF{t}.TrackObj.pf.z = tempDataList;
            TrackListPF{t}.TrackObj.pf = TrackListPF{t}.TrackObj.PredictMulti(TrackListPF{t}.TrackObj.pf);
            ValidationMatrix(t,:) = TrackListPF{t}.TrackObj.pf.Validation_matrix;
            tot_gate_area = tot_gate_area + TrackListPF{t}.TrackObj.pf.V_k;
        end
        % Compute New Track/False Alarm density
        bettaNTFA = sum(ValidationMatrix(:))/tot_gate_area;
        [TrackListPF] = JPDAF_EHM_PF_Update(TrackListPF, tempDataList, ValidationMatrix', bettaNTFA, 0);
        TrackNum = size(TrackList,2);
        %store Logs
        for j=1:TrackNum,
            Logs{j}.xV_ekf(:,i) = TrackListPF{j}.TrackObj.pf.xhk;
            st = [x1(i,j); y1(i,j)];
            Logs{j}.sV_ekf(:,i)= st;
            % Compute squared error
            Logs{j}.eV_ekf(:,i) = (TrackListPF{j}.TrackObj.pf.xhk(1:2,1) - st).*(TrackListPF{j}.TrackObj.pf.xhk(1:2,1) - st);
        end
        
        if (ShowPlots)
            if(i==1 || rem(i,SkipFrames+1)==0)
                % Plot data
                clf;
                 % Flip the image upside down before showing it
                imagesc([min_x max_x], [min_y max_y], flipud(img));

                % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

                hold on;
                for j=1:TrackNum,
                    h2 = plot(Logs{j}.sV_ekf(1,1:i),Logs{j}.sV_ekf(2,1:i),'b.-','LineWidth',1);
                    if j==2
                        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    end
                    h2 = plot(Logs{j}.sV_ekf(1,i),Logs{j}.sV_ekf(2,i),'bo','MarkerSize', 10);
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
                end
                h2 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
                for j=1:TrackNum,
                    colour = 'r';
                    if(j==2)
                       colour = 'c';
                    elseif (j==3)
                       colour = 'm';
                    end
                    h4 = plot(Logs{j}.xV_ekf(1,:),Logs{j}.xV_ekf(2,:),strcat(colour,'.-'),'LineWidth',1);
                    %h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
                    c_mean = mean(TrackListPF{j}.TrackObj.pf.particles,2);
                    c_cov = [std(TrackListPF{j}.TrackObj.pf.particles(1,:),TrackListPF{j}.TrackObj.pf.w')^2,0;0,std(TrackListPF{j}.TrackObj.pf.particles(2,:),TrackListPF{j}.TrackObj.pf.w')^2];
                    h2=plot_gaussian_ellipsoid(c_mean(1:2), c_cov);
                    set(h2,'color',colour);
                    set(h2,'LineWidth',1);
                    %plot(TrackListPF{j}.TrackObj.pf.particles(1,:),TrackListPF{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                    % set the y-axis back to normal.
                set(gca,'ydir','normal');
                str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
                title(str)
                xlabel('X position (m)')
                ylabel('Y position (m)')
                h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
                set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
                axis([0 10 0 10])


    %             subplot(1,3,2)
    %             hold on;
    %             set(gca, 'color', 'y')
    %             h1 = plot(1:N, Logs{1}.sV_ekf(1,:), 'b.-','LineWidth', 1);
    %             %set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    %             plot(1:N, Logs{2}.sV_ekf(1,:), 'b.-','LineWidth', 1);
    %             plot(1:N, Logs{1}.xV_ekf(1,:), 'r.-','LineWidth', 1);
    %             plot(1:N, Logs{2}.xV_ekf(1,:), 'c.-','LineWidth', 1);
    %             str = sprintf('Estimated state x_{1,k} vs. Time');
    %             title(str)
    %             xlabel('Time (s)')
    %             ylabel('X position (m)')
    %             h_legend = legend('Real', 'Target 1', 'Target 2');
    %             set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north', 'color', 'w');
    %             axis([0 180 2 10])
    %             
    %             subplot(1,3,3)
    %             hold on;
    %             hold on;
    %             set(gca, 'color', 'y')
    %             h1 = plot(1:N, Logs{1}.sV_ekf(2,:), 'b.-','LineWidth', 1);
    %             %set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    %             plot(1:N, Logs{2}.sV_ekf(2,:), 'b.-','LineWidth', 1);
    %             plot(1:N, Logs{1}.xV_ekf(2,:), 'r.-','LineWidth', 1);
    %             plot(1:N, Logs{2}.xV_ekf(2,:), 'c.-','LineWidth', 1);
    %             str = sprintf('Estimated state x_{1,k} vs. Time');
    %             title(str)
    %             xlabel('Time (s)')
    %             ylabel('Y position (m)')
    %             h_legend = legend('Real', 'Target 1', 'Target 2');
    %             set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north', 'color', 'w');
    %             axis([0 180 0 8])
                 pause(0.01)
            end
        end
    end
%     figure
%     subplot(2,1,1)
%     plot(1:N, Logs{1}.xV_ekf(1,:), 1:N, Logs{2}.xV_ekf(1,:));
%     str = sprintf('Estimated state x_{1,k} vs. Time');
%     title(str)
%     ylabel('Pendulum angle (rad)')
%     h_legend = legend('Meas', 'Real', 'EKF', 'UKF', 'PF');
%     set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'south');
%     
%     axis([0,5,-4,2])
%     subplot(2,1,2)
%     plot(1:N, Logs{1}.xV_ekf(2,:), 1:N, Logs{2}.xV_ekf(2,:));
%     str = sprintf('Estimated state x_{1,k} vs. Time');
%     title(str)
%     ylabel('Pendulum angle (rad)')
%     h_legend = legend('Meas', 'Real', 'EKF', 'UKF', 'PF');
%     set(h_legend,
    % Compute & Print RMSE
    RMSE = zeros(2, TrackNum);
    for i=1:TrackNum
%        RMSE(:,i) = sqrt(sum(Logs{i}.eV_ekf,2)/N)
    end
    RMSE_ekf = RMSE_ekf + RMSE;
end
RMSE_ekf = RMSE_ekf/SimNum
% figure
% hold on;
% h2 = plot(sV_ekf(1,1:N),sV_ekf(2,1:N),'b.-','LineWidth',1);
% h4 = plot(xV_ekf(1,1:N),xV_ekf(2,1:N),'r.-','LineWidth',1);
