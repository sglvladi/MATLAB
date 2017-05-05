% Number of Target Tracks
TrackNum = 2;

% Generate observations
%DataList = gen_obs_cluttered_multi(TrackNum, x_true, y_true);
RMSE_ekf = zeros(2, TrackNum);

% Number of simulations
SimNum = 1

% Show plots of data
ShowPlots = 0;
for i = 1:SimNum

    % Initiate KF parameters
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

    % Two-point difference initiation
    %for i = 1:TrackNum
    %    s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); DataList(1,i,2)-DataList(1,i,1); DataList(2,i,2)-DataList(2,i,1)]; %initial state
    %end
    %s.P_init = [q^2, 0, q^2, 0;
    %            0, q^2, 0, q^2;
    %            q^2, 0, 2*q^2, 0;
    %            0, q^2, 0, 2*q^2];                               % initial state covraiance

    % Single-point initiation
    Vmax = 0.4; % Max velocity = 0.4 m/s
    for i = 1:TrackNum
        s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); 0; 0]; %initial state
    end
    s.P_init = diag([q^2, q^2, (Vmax^2/3), (Vmax^2/3)]);


    % Process and Observation handles and covariances
    TrackObj.sys          = s.sys;
    TrackObj.obs          = s.obs;
    TrackObj.Q          = s.Q;
    TrackObj.R          = s.R;

    % initial covariMeasInd, TrackIndance assumed same for all tracks
    TrackObj.P          = s.P_init;

    N=size(DataList,3) ;                                    % total dynamic steps

    % Instantiate EKF and vars to store output
    Logs = [];
    Log.xV_ekf = zeros(n,N);          %estmate        % allocate memory
    Log.PV_ekf = zeros(1,N);
    %PV_ekf = cell(1,N);     % use to display ellipses 
    Log.sV_ekf = zeros(n/2,N);          %actual
    Log.zV_ekf = zeros(n/2,N);
    Log.eV_ekf = zeros(n/2,N);

    % Initiate Tracklist
    TrackList = [];
    for i=1:TrackNum,

        TrackObj.x          = s.x_init(:,i)
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
    for i=2:N

        i
        if(i>200)
            i
        end
        [TrackList, ValidationMatrix, bettaNTFA] = Observation_Association(TrackList, DataList(:,:,i), ekf);
        %TrackList = JPDAF_EKF_Update(TrackList, DataList(:,:,i), ValidationMatrix', bettaNTFA);
        TrackList = JPDAF_UKF_Update(TrackList, DataList(:,:,i), ValidationMatrix', bettaNTFA);
        %store Logs
        for j=1:TrackNum,
            Logs{j}.xV_ekf(:,i) = TrackList{j}.TrackObj.x;
            st = [x_true(i,j); y_true(i,j)];
            Logs{j}.sV_ekf(:,i)= st;
            % Compute squared error
            Logs{j}.eV_ekf(:,i) = (TrackList{j}.TrackObj.x(1:2,1) - st).*(TrackList{j}.TrackObj.x(1:2,1) - st);
        end
        
        if (ShowPlots)
            % Plot data
            clf; % Flip the image upside down before showing it
            imagesc([min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            for j=1:TrackNum,
                h2 = plot(Logs{j}.sV_ekf(1,1:i),Logs{j}.sV_ekf(2,1:i),'b.-','LineWidth',1);
                h2 = plot(Logs{j}.sV_ekf(1,i),Logs{j}.sV_ekf(2,i),'bo','MarkerSize', 20);
            end
            h2 = plot(DataList(1,:,i),DataList(2,:,i),'k*','MarkerSize', 20);
            for j=1:TrackNum,
                colour = 'r';
                if(j==2)
                    colour = 'c';
                elseif (j==3)
                    colour = 'm';
                end
                h4 = plot(Logs{j}.xV_ekf(1,1:i),Logs{j}.xV_ekf(2,1:i),strcat(colour,'.-'),'LineWidth',1);
                h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 20);
            end
                % set the y-axis back to normal.
            set(gca,'ydir','normal');
            axis([0 10 0 10])
            pause(0.1)
        end
    end

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
