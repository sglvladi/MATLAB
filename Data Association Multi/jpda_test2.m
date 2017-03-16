% Number of Target Tracks
TrackNum = 3;

% Generate observations
[DataList,x1,y1] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, 0.5, 2, 10,1);
RMSE_ekf = zeros(2, TrackNum);

% Number of simulations
SimNum = 1

% Show plots of data
ShowPlots = 1;
for sim = 1:SimNum

    % Initiate KF parameters
     n=4;      %number of state
     q=0.01;    %std of process 
     r=0.25;    %std of measurement
     s.Q=[1^3/3, 0, 1^2/2, 0;  0, 1^3/3, 0, 1^2/2; 1^2/2, 0, 1, 0; 0, 1^2/2, 0, 1]*10*q^2; % covariance of process
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
        s.x_init(:,i)=DataList{1}(1,i); DataList{1}(2,i); 0; 0; %initial state
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
% 
%     % Initiate Tracklist
     TrackList = [];
    for i=1:TrackNum,

        TrackObj.x          = s.x_init(:,i);
        TrackList{i}.TrackObj = TrackObj;
        Logs{i}             = Log;
        Logs{i}.xV_ekf(:,1) = s.x_init(:,i);
        Logs{i}.sV_ekf(:,1) = st(:,i);
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
        [TrackList, ValidationMatrix, bettaNTFA] = Observation_Association(TrackList, tempDataList, ekf);
        %TrackList = JPDAF_EKF_Update(TrackList, DataList(:,:,i), ValidationMatrix', bettaNTFA);
        [TrackList, betta] = JPDAF_EHM_Update(TrackList, tempDataList, ValidationMatrix', bettaNTFA);
        %TrackList = Track_InitConfDel(TrackList,tempDataList,ValidationMatrix',bettaNTFA, betta);
        
        
        TrackNum = size(TrackList,2);
        %store Logs
        for j=1:TrackNum,
            Logs{j}.xV_ekf(:,i) = TrackList{j}.TrackObj.x;
            %st = [x1(i,j); y1(i,j)];
           % Logs{j}.sV_ekf(:,i)= st;
            % Compute squared error
            %Logs{j}.eV_ekf(:,i) = (TrackList{j}.TrackObj.x(1:2,1) - st).*(TrackList{j}.TrackObj.x(1:2,1) - st);
        end
        
        if (ShowPlots)
            % Plot data
            clf;
             % Flip the image upside down before showing it
            imagesc([min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
%             for j=1:TrackNum,
%                 %h2 = plot(Logs{j}.sV_ekf(1,1:i),Logs{j}.sV_ekf(2,1:i),'b.-','LineWidth',1);
%                 if j==2
%                     set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%                 end
%                 %h2 = plot(Logs{j}.sV_ekf(1,i),Logs{j}.sV_ekf(2,i),'bo','MarkerSize', 10);
%                 set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
%             end
            h2 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
            for j=1:TrackNum,
                colour = 'r';
                if(j==2)
                   colour = 'c';
                elseif (j==3)
                   colour = 'm';
                end
                h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'.-'),'LineWidth',1);
                h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
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
             pause(0.1)
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
        RMSE(:,i) = sqrt(sum(Logs{i}.eV_ekf,2)/N)
    end
    RMSE_ekf = RMSE_ekf + RMSE;
end
RMSE_ekf = RMSE_ekf/SimNum
% figure
% hold on;
% h2 = plot(sV_ekf(1,1:N),sV_ekf(2,1:N),'b.-','LineWidth',1);
% h4 = plot(xV_ekf(1,1:N),xV_ekf(2,1:N),'r.-','LineWidth',1);
