%   jpda_test.m                                  Author: Lyudmil Vladimirov
%   ======================================================================>
%   Functionality: Main loop to perform JPDA testing
%   
%   Input: x_true, y_true - Matrices (m x n) of x and y coordinates 
%                           (m: time, n: targets) 
%   
%   Dependencies: gen_obs_cluttered_multi.m, UKalmanFilter.m,
%                 Observation Association.m, JPDAF_UKF_Update.m 
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Number of targets to be tracked (e.g. n)
TrackNum = size(x_true, 2);

% Generate measurements (including clutter) from ground truth
DataList = gen_obs_cluttered_multi(TrackNum, x_true, y_true); 

% Create variable to store RMSE for each track
RMSE_ukf = zeros(2, TrackNum);

% Number of simulation runs
SimNum = 1

% Plotting variables
ShowPlots = 1; % Enable/Disable plots
SkipFrames = 5; % Skip n frames between plots

% for each simulation run
for sim = 1:SimNum

    %% Initiate KF parameters
    n=4;      %number of state
    q=0.01;    %std of process 
    r=0.25;    %std of measurement
    s.Q=[1^3/3, 0, 1^2/2, 0;  0, 1^3/3, 0, 1^2/2; 1^2/2, 0, 1, 0; 0, 1^2/2, 0, 1]*10*q^2; % covariance of process
    s.R=r^2*eye(n/2);        % covariance of measurement  
    s.sys=(@(x)[x(1)+ x(3); x(2)+x(4); x(3); x(4)]);  % assuming measurements arrive 1 per sec
    s.obs=@(x)[x(1);x(2)];                               % measurement equation
    st=[x_true(1,:);y_true(1,:)];                                % initial state
    s.x_init = [];
    s.P_init = [];

    %% Perform Track Initiations
    
    % 1) Two-point difference initiation
    %for i = 1:TrackNum
    %    s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); DataList(1,i,2)-DataList(1,i,1); DataList(2,i,2)-DataList(2,i,1)]; %initial state
    %end
    %s.P_init = [q^2, 0, q^2, 0;
    %            0, q^2, 0, q^2;
    %            q^2, 0, 2*q^2, 0;
    %            0, q^2, 0, 2*q^2];     % initial state covraiance

    % 2)Single-point initiation
    Vmax = 0.4; % Max velocity = 0.4 m/s
    for i = 1:TrackNum
        s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); 0; 0]; %initial state
    end
    s.P_init = diag([q^2, q^2, (Vmax^2/3), (Vmax^2/3)]);

    N=size(DataList,3) ;   % total dynamic steps

    %% Instantiate vars to store output
    Logs = [];
    Log.xV_ukf = zeros(n,N);    % mean estimates
    Log.PV_ukf = zeros(n,n,N);    % covariance estimates    
    Log.sV_ukf = zeros(n/2,N);  % ground truth
    Log.zV_ukf = zeros(n/2,N);  % measurements
    Log.eV_ukf = zeros(n/2,N);  % RMSE
    
    %% Process and Observation handles and covariances
    TrackObj.sys          = s.sys;
    TrackObj.obs          = s.obs;
    TrackObj.Q          = s.Q;
    TrackObj.R          = s.R;
    TrackObj.P          = s.P_init;

    %% Initiate Tracklist and Log for each target
    TrackList = [];
    for i=1:TrackNum,
        
        TrackObj.x          = s.x_init(:,i);
        TrackList{i}.TrackObj = TrackObj;
        
        Logs{i}             = Log;
        Logs{i}.xV_ukf(:,1) = s.x_init(:,i);
        Logs{i}.PV_ukf(:,:,1) = s.P_init;
        Logs{i}.sV_ukf(:,1) = st(:,i);
    end;

    %% Create UKF instance to perform Track prediction in Observation Association function 
    ukf = UKalmanFilter(TrackObj, 0.5, 1, 2);

    %% Prepare background image from plot
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
    
    %% Simulation
    for i=1:N

        fprintf('Iteration = %d/%d\n',i,N);
        %tempDataList = DataList(:,:,i);
        %tempDataList( :, ~any(tempDataList,1) ) = []; 
        
        % Run Observation Association to perform track update for all targets
        %   and generate the Validation matrix.
        [TrackList, ValidationMatrix, bettaNTFA] = Observation_Association(TrackList, DataList(:,:,i), ukf);
        
        % JPDAF_UKF_Update provides JPDAF updates for all targets, using a 
        %   UKF.
        TrackList = JPDAF_EHM_Update(TrackList, DataList(:,:,i), ValidationMatrix', bettaNTFA);
        %TrackList = Track_InitConfDel(TrackList,tempDataList,ValidationMatrix',bettaNTFA, betta);
        
        % Update Logs
        for j=1:TrackNum
            Logs{j}.xV_ukf(:,i) = TrackList{j}.TrackObj.x;
            Logs{j}.PV_ukf(:,:,i) = TrackList{j}.TrackObj.P;
            st = [x_true(i,j); y_true(i,j)];
            Logs{j}.sV_ukf(:,i)= st;
            % Compute squared error
            Logs{j}.eV_ukf(:,i) = (TrackList{j}.TrackObj.x(1:2,1) - st).*(TrackList{j}.TrackObj.x(1:2,1) - st);
        end
        
        % Visualise the process (Set ShowPlots=0 to disable)
        if (ShowPlots)
            % Plot data
            
            if(i==1 || rem(i,SkipFrames+1)==0)
                clf;
                 
                 % Flip the image upside down before showing it
                imagesc([min_x max_x], [min_y max_y], flipud(img));
                hold on;
                
                % Draw ground truth for all targets
                for j=1:TrackNum,
                    h2 = plot(Logs{j}.sV_ukf(1,1:i),Logs{j}.sV_ukf(2,1:i),'b.-','LineWidth',1);
                    if j~=1
                        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    end
                    h2 = plot(Logs{j}.sV_ukf(1,i),Logs{j}.sV_ukf(2,i),'bo','MarkerSize', 10);
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
                end
                
                % Draw all measurements
                h2 = plot(DataList(1,:,i),DataList(2,:,i),'k*','MarkerSize', 10);
                
                % Draw estimated position mean and covariance
                for j=1:TrackNum,
                    colour = 'r';
                    if(j==2)
                        colour = 'c';
                    elseif (j==3)
                        colour = 'm';
                    end
                    h4 = plot(Logs{j}.xV_ukf(1,1:i),Logs{j}.xV_ukf(2,1:i),strcat(colour,'.-'),'LineWidth',1);
                    %h4 = plot(Logs{j}.xV_ukf(1,i),Logs{j}.xV_ukf(2,i),strcat(colour,'o'),'MarkerSize', 10);
                    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    h2=plot_gaussian_ellipsoid(Logs{j}.xV_ukf(1:2,i), Logs{j}.PV_ukf(1:2,1:2,i));
                    set(h2,'color',colour);
                    set(h2,'LineWidth',1);
                end
                
                % set the y-axis back to normal.
                set(gca,'ydir','normal');
                str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
                title(str);
                xlabel('X position (m)');
                ylabel('Y position (m)');
                
                if TrackNum==1
                    h_legend = legend('Real', 'Meas', 'Target 1');
                elseif TrackNum == 2
                    h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
                elseif TrackNum == 3
                    h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2', 'Target 3');
                end
                
                set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
                axis([0 10 0 10]);
                pause(0.01)
            end
        end
    end

    % Compute & Print RMSE
    RMSE = zeros(2, TrackNum);
    for i=1:TrackNum
        RMSE(:,i) = sqrt(sum(Logs{i}.eV_ukf,2)/N);
    end
    RMSE_ukf = RMSE_ukf + RMSE;
end
RMSE_ukf = RMSE_ukf/SimNum

%% Draw history/end plot
clf;
% Flip the image upside down before showing it
imagesc([min_x max_x], [min_y max_y], flipud(img));

hold on;
 for j=1:TrackNum,
    h2 = plot(Logs{j}.sV_ukf(1,:),Logs{j}.sV_ukf(2,:),'b.-','LineWidth',1);
    if j>=2
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    h2 = plot(Logs{j}.sV_ukf(1,end),Logs{j}.sV_ukf(2,end),'bo','MarkerSize', 10);
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
end
h2 = plot(DataList(1,:,end),DataList(2,:,end),'k*','MarkerSize', 10);
for j=1:TrackNum,
    colour = 'r';
    if(j==2)
        colour = 'c';
    elseif (j==3)
        colour = 'm';
    end
    h4 = plot(Logs{j}.xV_ukf(1,:),Logs{j}.xV_ukf(2,:),strcat(colour,'.-'),'LineWidth',1);
    %h4 = plot(Logs{j}.xV_ukf(1,end),Logs{j}.xV_ukf(2,end),strcat(colour,'o'),'MarkerSize', 10);
    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h2=plot_gaussian_ellipsoid(Logs{j}.xV_ukf(1:2,end), Logs{j}.PV_ukf(1:2,1:2,end));
    set(h2,'color',colour);
    set(h2,'LineWidth',1);
end
    % set the y-axis back to normal.
set(gca,'ydir','normal');
str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
title(str);
xlabel('X position (m)');
ylabel('Y position (m)');
if TrackNum==1
    h_legend = legend('Real', 'Meas', 'Target 1');
elseif TrackNum == 2
    h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
elseif TrackNum == 3
    h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2', 'Target 3');
end
set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
axis([0 10 0 10]);
