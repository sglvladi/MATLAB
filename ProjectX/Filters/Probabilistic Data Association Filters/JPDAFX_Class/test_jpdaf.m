%% Plot settings
ShowPlots = 1;
ShowArena = 1;
SkipFrames = 0;

%% Constant Velocity model
config_cv.q = 0.0001;
CVmodel = ConstantVelocityModelX(config_cv);

%% Constant Heading Model
config_ch.q_vel = 0.01;
config_ch.q_head = 0.3;
CHmodel = ConstantHeadingModelX(config_ch);

%% 2D linear-Gaussian Observation Model
config_pom.r = 0.3;
POmodel = PositionalObsModel(config_pom);

%% Assign PF parameter values
config_pf.k               = 1;                   % initial iteration number
config_pf.Np              = 5000;                 % number of particles
config_pf.resampling_strategy = 'systematic_resampling';

%% Initiate Kalman Filter
config_kf.k = 1;


%% Set TrackNum
TrackNum = 3;

%% Generate DataList
[DataList,x1,y1] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, 0.3, 2, 50, 1);

%% Get GroundTruth
for i=1:TrackNum
    GroundTruth{i} = [x_true(:,i), y_true(:,i)]; % ith target's GroundTruth
end

%% Initiate TrackList
for i=1:TrackNum
    config_kf.x = [GroundTruth{i}(1,1);GroundTruth{i}(1,2); 0; 0];
    config_kf.P = CVmodel.config.Q(1);
    %FilterList{1}.Filter = KalmanFilterX(config_kf, CVmodel, obs_model);
    config_pf.gen_x0 =  @(Np) [mvnrnd(repmat([GroundTruth{i}(1,1),GroundTruth{i}(1,2)],Np,1),POmodel.config.R(1)), rand(Np,1), 2*pi*rand(Np,1)];
    TrackList{i}.TrackObj = EKalmanFilterX(config_kf, CHmodel, POmodel);%ParticleFilterX(config_pf, CHmodel, POmodel);  
end

%% Initiate PDAF parameters
config_jpdaf.Y = DataList{1}(:,:);  
config_jpdaf.TrackList = TrackList;
config_jpdaf.PD = 0.8;
config_jpdaf.PG = 0.998; %0.999;
config_jpdaf.GateLevel = 10;%14;

%% Instantiate JPDAF
jpdaf = JPDAFX(config_jpdaf);

%% Instantiate Log to store output
N=size(DataList,2);
Logs = cell(1, 5); % 4 tracks
N = size(x_true,1)-2;
for i=1:TrackNum
    Logs{i}.sv = zeros(4,N);
    Logs{i}.xV = zeros(4,N);          %estmate        % allocate memory
    Logs{i}.err = zeros(2,N);
    Logs{i}.pos_err = zeros(1,N);
    Logs{i}.exec_time = 0;
    Logs{i}.filtered_estimates = cell(1,N);
end

% Create figure windows
if(ShowPlots)
    if(ShowArena)
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
    end
    
    figure('units','normalized','outerposition',[0 0 .5 1])
    ax(1) = gca;
end

exec_time = 0;
for i = 1:N
    % Remove null measurements   
    jpdaf.config.k = 1;
    jpdaf.config.Y = DataList{i}(:,:);
    jpdaf.config.Y( :, ~any(jpdaf.config.Y,1) ) = [];
    
    for j=1:TrackNum
        %Logs{j}.xV_ekf(:,i) = jpdaf.config.TrackList{j}.TrackObj.config.x;
        st = [x1(i,j); y1(i,j)];
        Logs{j}.sV_ekf(:,i)= st;
        %Logs{j}.eV_ekf(:,i) = (jpdaf.config.TrackList{j}.TrackObj.config.x(1:2,1) - st).*(jpdaf.config.TrackList{j}.TrackObj.config.x(1:2,1) - st);
    end
    tic;
    jpdaf.Predict();
    %exec_time = 
    if (ShowPlots && i>1)
        if(i==1 || rem(i,SkipFrames+1)==0)
            % Plot data
            clf;
             % Flip the image upside down before showing it
            imagesc([min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            for j=1:TrackNum
                h2 = plot(Logs{j}.sV_ekf(1,1:i),Logs{j}.sV_ekf(2,1:i),'b.-','LineWidth',1);
                if j==2
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                h2 = plot(Logs{j}.sV_ekf(1,i),Logs{j}.sV_ekf(2,i),'bo','MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            end
            h2 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
            for j=1:TrackNum
                colour = 'r';
                if(j==2)
                   colour = 'c';
                elseif (j==3)
                   colour = 'm';
                end
                h4 = plot(Logs{j}.xV_ekf(1,:),Logs{j}.xV_ekf(2,:),strcat(colour,'.-'),'LineWidth',1);
                %h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
                if(isa(jpdaf.config.TrackList{j}.TrackObj,'ParticleFilterX'))
                    c_mean = sum(repmat(jpdaf.config.TrackList{j}.TrackObj.config.w,size(jpdaf.config.TrackList{j}.TrackObj.config.particles,1),1).*jpdaf.config.TrackList{j}.TrackObj.config.particles,2);
                    c_cov = [std(jpdaf.config.TrackList{j}.TrackObj.config.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.config.w)^2,0;0,std(jpdaf.config.TrackList{j}.TrackObj.config.particles(2,:),jpdaf.config.TrackList{j}.TrackObj.config.w)^2];
                else
                    c_mean = jpdaf.config.TrackList{j}.TrackObj.config.x_pred;
                    c_cov = jpdaf.config.TrackList{j}.TrackObj.config.P_pred(1:2,1:2);
                end
                h2=plot_gaussian_ellipsoid(c_mean(1:2), c_cov);
                set(h2,'color',colour);
                set(h2,'LineWidth',1);
                %plot(jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
                % set the y-axis back to normal.
            set(gca,'ydir','normal');
            str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
            title(str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
%            h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
%            set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
            axis([0 10 0 10])
            pause(0.01)
        end
    end
    
    jpdaf.Update();
        
    %store Logs
    for j=1:TrackNum
        Logs{j}.xV_ekf(:,i) = jpdaf.config.TrackList{j}.TrackObj.config.x;
        st = [x1(i,j); y1(i,j)];
        Logs{j}.sV_ekf(:,i)= st;
        Logs{j}.eV_ekf(:,i) = (jpdaf.config.TrackList{j}.TrackObj.config.x(1:2,1) - st).*(jpdaf.config.TrackList{j}.TrackObj.config.x(1:2,1) - st);
    end

    if (ShowPlots)
        if(i==1 || rem(i,SkipFrames+1)==0)
            % Plot data
            clf;
             % Flip the image upside down before showing it
            imagesc([min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            for j=1:TrackNum
                h2 = plot(Logs{j}.sV_ekf(1,1:i),Logs{j}.sV_ekf(2,1:i),'b.-','LineWidth',1);
                if j==2
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                h2 = plot(Logs{j}.sV_ekf(1,i),Logs{j}.sV_ekf(2,i),'bo','MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            end
            h2 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
            for j=1:TrackNum
                colour = 'r';
                if(j==2)
                   colour = 'c';
                elseif (j==3)
                   colour = 'm';
                end
                h4 = plot(Logs{j}.xV_ekf(1,:),Logs{j}.xV_ekf(2,:),strcat(colour,'.-'),'LineWidth',1);
                %h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
                if(isa(jpdaf.config.TrackList{j}.TrackObj,'ParticleFilterX'))
                    c_mean = sum(repmat(jpdaf.config.TrackList{j}.TrackObj.config.w,size(jpdaf.config.TrackList{j}.TrackObj.config.particles,1),1).*jpdaf.config.TrackList{j}.TrackObj.config.particles,2);
                    c_cov = [std(jpdaf.config.TrackList{j}.TrackObj.config.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.config.w)^2,0;0,std(jpdaf.config.TrackList{j}.TrackObj.config.particles(2,:),jpdaf.config.TrackList{j}.TrackObj.config.w)^2];
                else
                    c_mean = jpdaf.config.TrackList{j}.TrackObj.config.x;
                    c_cov = jpdaf.config.TrackList{j}.TrackObj.config.P(1:2,1:2);
                end
                h2=plot_gaussian_ellipsoid(c_mean(1:2), c_cov);
                set(h2,'color',colour);
                set(h2,'LineWidth',1);
                %plot(jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
                % set the y-axis back to normal.
            set(gca,'ydir','normal');
            str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
            title(str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
%            h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
%            set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
            axis([0 10 0 10])
            pause(0.01)
        end
    end
end