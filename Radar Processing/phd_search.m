%% Plot settings
ShowPlots = 1;
ShowPredict = 0;
ShowData = 1;
SkipFrames = 0;

%% Recording Settings
Record = 1;
clear F
clear M

%% Initiate PF parameters
nx = 4;      % number of state dims
nu = 4;      % size of the vector of process noise
nv = 2;      % size of the vector of observation noise
q  = 1;   % process noise density (std)
r  = 20;    % observation noise density (std)
lambdaV = 10; % mean number of clutter points 
V_bounds = [-2000 -800 2000 3000]; %[-2.5 .2 -3 3]; [-2 -.800 2 3] [-.700 -.400 -.700 .400]; % [x_min x_max y_min y_max]
V = (abs(V_bounds(2)-V_bounds(1))*abs(V_bounds(4)-V_bounds(3)));
% Prior PDF generator
gen_x0_cch = @(Np) mvnrnd(repmat([0,0,0,0],Np,1),diag([q^2, q^2, 100, 100]));
% Process equation x[k] = sys(k, x[k-1], u[k]);
sys_cch = @(k, xkm1, uk) [xkm1(1,:)+k*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+k*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ k*uk(:,3)'; xkm1(4,:) + k*uk(:,4)'];
% PDF of process noise generator function
gen_sys_noise_cch = @(u) mvnrnd(zeros(size(u,2), nu), diag([0,0,q^2,0.16^2])); 
% Observation equation y[k] = obs(k, x[k], v[k]);
obs = @(k, xk, vk) [xk(1)+vk(1); xk(2)+vk(2)];                  % (returns column vector)
% PDF of observation noise and noise generator function
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
pf.particles       = zeros(4, pf.Np); % particles
pf.resampling_strategy = 'systematic_resampling';
pf.sys = sys_cch;
pf.particles = zeros(nx, pf.Np); % particles
pf.gen_x0 = gen_x0_cch(pf.Np);
pf.obs = p_yk_given_xk;
pf.obs_model = @(xk) [xk(1,:); xk(2,:)];
pf.R = cov_v;
pf.clutter_flag = 1;
pf.multi_flag = 1;
pf.sys_noise = gen_sys_noise_cch;
pf.V_k = 0;

%% Set TrackNum
TrackNum = 0;
TrueTracks = 3;

%% Generate DataList                       (TrackNum, x_true, y_true, R, R_clutter, lambdaV, Iter)                   
%[DataList,x1,y1] = gen_obs_cluttered_multi2(TrueTracks, x_true, y_true, r, 2, lambdaV, 1);

%% Get GroundTruth
% for i=1:TrueTracks
%     GroundTruth{i} = [x_true(:,i), y_true(:,i)]; % ith target's GroundTruth
% end

%% Initiate TrackList
% % for i=1:TrackNum,
%      pf.gen_x0 = @(Np) [mvnrnd(repmat([-538.400000000000,-144],Np,1),cov_v), 5^2*rand(Np,1), 2*pi*rand(Np,1)];
%      pf.ExistProb = 0.8;
%      TrackList{1}.TrackObj = ParticleFilterMin2(pf);
%      pf.gen_x0 = @(Np) [mvnrnd(repmat([-558.400000000000, 97.4000000000000],Np,1),cov_v), 5^2*rand(Np,1), 2*pi*rand(Np,1)];
%      pf.ExistProb = 0.8;
%      TrackList{2}.TrackObj = ParticleFilterMin2(pf);
% end;

%% Initiate PDAF parameters
Par = [];
Par.Filter = ParticleFilterMin2(pf);
Par.DataList = DataList{1}(:,:);
%Par.GroundTruth = GroundTruth;
Par.TrackList = [];
Par.PD = 0.5;
Par.PG = 0.998;
Par.GateLevel = 15;
Par.Pbirth = 0.01;
Par.Pdeath = 0.005;
Par.SimIter = 1000;

%% Assign PHD parameter values
par.k               = 1;                                                    % initial iteration number
par.Np              = 100000;                                                % number of particles
par.resampling_strategy = 'systematic_resampling';                          % resampling strategy
par.birth_strategy = 'mixture';                                           %  
par.sys = @(k, xkm1, uk) [xkm1(1,:)+k*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+k*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ k*uk(:,3)'; xkm1(4,:) + k*uk(:,4)']; % CH model
%par.gen_x0 = @(Np)[(2000+200)*rand(Np,1)-2000,6000*rand(Np,1)-3000, mvnrnd(zeros(Np,1), 2^2), 2*pi*rand(Np,1)]; % Uniform position and heading, Gaussian speed
par.gen_x0 = @(Np)[abs(V_bounds(2)-V_bounds(1))*rand(Np,1)+V_bounds(1),abs(V_bounds(4)-V_bounds(3))*rand(Np,1)+V_bounds(3), 2^2*rand(Np,1), 2*pi*rand(Np,1)]; % Uniform position and heading, Gaussian speed
par.particles = par.gen_x0(par.Np)';                                        % Generate inital particles as per gen_x0
par.gen_x1 = @(obs_mean, Np) [mvnrnd(repmat(obs_mean,1,Np)', cov_v), mvnrnd(zeros(Np,1), 2^2), 2*pi*rand(Np,1)];
%par.particles = par.gen_x0([0;0], par.Np)'; % Generate inital particles as per gen_x0
par.w = repmat(1/par.Np, par.Np, 1)';                                       % Uniform weights
par.likelihood = @(k, yk, xk) mvnpdf(yk, xk, cov_v);                        % Likelihood model p(y|x)
par.obs_model = @(xk) [xk(1,:); xk(2,:)];                                   % Observation model (no noise)
par.sys_noise = gen_sys_noise_cch;                                          % System noise
par.Pbirth = 0.01;                                                           % Birth probability
par.Pdeath = 0.005;                                                         % Death probability
par.J_k = 1000;                                                             % Number of birth particles (births by "expansion")
par.PD = 0.5;                                                               % Probability of detection
par.lambda = lambdaV/V;                                                % Mean clutter per unit area
par.Np_conf = pf.Np;                                                        % Number of particles for confirmed tracks
par.P_conf = 0.90;                                                           % Confirmation probability
par.type = 'search';                                                        % Search PHD filter
par.R = cov_v;

myphd = SMC_PHD(par);

%% Instantiate JPDAF
jpdaf = JPDAF(Par);
% jpdaf.config.TrackList = TrackList;

%% Instantiate Log to store output
N=size(DataList,2);
Logs = [];
Log.xV_ekf = zeros(nx,N);          %estmate        % allocate memory
Log.PV_ekf = zeros(1,N);
Log.sV_ekf = zeros(nx/2,N);          %actual
Log.zV_ekf = zeros(nx/2,N);
Log.eV_ekf = zeros(nx/2,N);


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

figure('units','normalized','outerposition',[0 0 .5 1])
ax(1) = gca;
figure('units','normalized','outerposition',[.5 0 .5 1])
ax(2) = gca;

exec_time = 0;
for i = 1:N
    tic;
    fprintf('\nIteraration: %d/%d\n', i, N);
    
    % Remove null measurements   
    DataList_k = DataList{i}(:,:);
    DataList_k( :, ~any(DataList_k,1) ) = [];
    
    % Change JPDAF parameters
    jpdaf.config.DataList = DataList_k; % New observations
    
    % Change PHD filter parameters
    myphd.config.lambda = size(DataList_k,2)/V;
    myphd.config.k = i; % Time index
    myphd.config.z = DataList_k; % New observations
    
    for t = 1:length(jpdaf.config.TrackList)
        if(i==1)
            myphd.config.k = i; % Time index
            jpdaf.config.TrackList{t}.TrackObj.pf.k = i;
        else
            myphd.config.k = Dt(i-1); % Time index
            jpdaf.config.TrackList{t}.TrackObj.pf.k = Dt(i-1);
        end
    end
    
    
    % 1) Predict the confirmed tracks
    jpdaf.Predict();
        
    % 2) Predict the PHD filter
    myphd.Predict();
    if(ShowPlots&&ShowPredict)
        cla(ax(1));
        if(ShowData)
            h2 = plot(ax(1), DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
        end
        for j=1:jpdaf.config.TrackNum
            colour = 'r';
            if(j==2)
               colour = 'c';
            elseif (j==3)
               colour = 'm';
            end
            h4 = plot(ax(1), jpdaf.config.TrackList{j}.TrackObj.pf.xhk(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.xhk(2,:),strcat(colour,'.-'),'LineWidth',1);
            %h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
            c_mean = mean(jpdaf.config.TrackList{j}.TrackObj.pf.particles,2);
            c_cov = [std(jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.w')^2,0;0,std(jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),jpdaf.config.TrackList{j}.TrackObj.pf.w')^2];
            h2=plot_gaussian_ellipsoid(c_mean(1:2), c_cov, 1, [], ax(1));
            set(h2,'color',colour);
            set(h2,'LineWidth',1);
            %plot(ax(1),jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
            set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            text(ax(1),c_mean(1)+.020,c_mean(2)-.05,int2str(j));
            text(ax(1),c_mean(1)+.020,c_mean(2)-.15,num2str(jpdaf.config.TrackList{j}.TrackObj.pf.ExistProb, 2));
        end
            % set the y-axis back to normal.
        %set(ax(1),'ydir','normal');
        str = sprintf('Visualisation of tracking process');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
    %            h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
    %            set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
        axis(ax(1),V_bounds)
        %pause(0.1)
    
        % Plot PHD
        cla(ax(2), 'reset');
        [bandwidth,density,X,Y]=kde2d(myphd.config.particles(1:2,:)');
        %contour3(X,Y,density,50);
        surf(ax(2),X,Y,density);
        hold on;
        plot(ax(2), myphd.config.particles(1,:), myphd.config.particles(2,:), '.')
        hold on;
        plot(ax(2), myphd.config.z(1,:), myphd.config.z(2,:), 'y*');
        axis(ax(2), V_bounds);
        str = sprintf('Visualisation of PHD search track density');
        xlabel(ax(2),'X position (m)')
        ylabel(ax(2),'Y position (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        %pause(0.01)
    end
    
    % 3) Compute rhi as given by Eq. (16) in [2] 
    rhi = zeros(1, size(DataList_k,2));
    for m = 1:size(DataList_k,2)
        rhi_tmp = 1;
        if(jpdaf.config.betta>-1) % Check if beta exists
            for t = 1:TrackNum
                rhi_tmp = rhi_tmp*(1-jpdaf.config.betta(t,m+1));
            end
        end
        rhi(m) = rhi_tmp;
    end
    myphd.config.rhi = rhi;
    
    % 4) Update the confirmed track
    jpdaf.Update();
    
    % 5) Update PHD filter
    myphd.Update();
    
    % 6) Initiate tracks
    for j = 1:length(myphd.config.NewTracks)
        pf.particles = myphd.config.NewTracks{j}.particles;
        pf.w = myphd.config.NewTracks{j}.w;
        pf.ExistProb = myphd.config.NewTracks{j}.ExistProb;
        jpdaf.config.TrackList{end+1}.TrackObj = ParticleFilterMin2(pf);
        jpdaf.config.TrackNum = jpdaf.config.TrackNum + 1;
    end
    % 7) Delete tracks
     del_tracks = 0;
     del_flag = 0;
    for t = 1:length(jpdaf.config.TrackList)
        if(jpdaf.config.TrackList{t}.TrackObj.pf.ExistProb<0.1 || jpdaf.config.TrackList{t}.TrackObj.pf.V_k > 60000)
            jpdaf.config.TrackList{t} = [];
            del_tracks = del_tracks + 1;
            del_flag = 1;
        end
    end
    if(del_flag)
        jpdaf.config.TrackList = jpdaf.config.TrackList(~cellfun('isempty',jpdaf.config.TrackList));
        jpdaf.config.TrackNum = jpdaf.config.TrackNum - del_tracks;
    end
    TrackNum = size(jpdaf.config.TrackList,2);
    
    exec_time = exec_time + toc;
%     if(TrackNum>0)
%         % Extract valid tracks
%         ValidTracks = 0;
%         TrackDists = zeros(TrackNum, TrueTracks);
%         for j=1:TrackNum
%             for t=1:TrueTracks
%                 if(sqrt((jpdaf.config.TrackList{j}.TrackObj.pf.xhk(1) - x1(i,t))^2 + (jpdaf.config.TrackList{j}.TrackObj.pf.xhk(2)-y1(i,t))^2)<3*r)
%                     TrackDists(j,t) = sqrt((jpdaf.config.TrackList{j}.TrackObj.pf.xhk(1) - x1(i,t))^2 + (jpdaf.config.TrackList{j}.TrackObj.pf.xhk(2)-y1(i,t))^2);
%                 else
%                     TrackDists(j,t) = 10^15 + t;
%                 end
%             end
%         end
%         [B,I] = sort(TrackDists,1);
%         sorted_inds = I;
%         for t=1:TrueTracks
%             if TrackDists(sorted_inds(1,t),t)>=10^15
%                 sorted_inds(1,t) = TrackNum+t;
%             end
%         end
%         if(length(unique(sorted_inds(1,:)))<length(sorted_inds(1,:)))
%             %for j=1:TrackNum
%     %             for t=1:TrueTracks
%     %                 min_inds_tmp = min_inds;
%     %                 min_inds_tmp(t) = [];
%     %                 if(~isempty(intersect(min_inds(t), min_inds_tmp)))
%             %error('Sort error');
%             [B,Ia,Ib] = intersect(unique(sorted_inds(1,:)), sorted_inds(1,:));
%             sorted_inds(1,Ib(1)) =  sorted_inds(2,Ib(1));
%         end
% 
%         for t=1:TrueTracks
%             if sorted_inds(1,t)<=TrackNum
%                 ValidTracks = ValidTracks+1;
%             end
%         end
% 
%         ValidTrackNum_log(i) = ValidTracks;
%     else
%         ValidTrackNum_log(i) = 0;
%     end
    
    %store Logs
    for j=1:jpdaf.config.TrackNum
        try
            Logs{j}.xV_ekf(:,end+1) = jpdaf.config.TrackList{j}.TrackObj.pf.xhk;
        catch
            Logs{j}.xV_ekf(:,1) = jpdaf.config.TrackList{j}.TrackObj.pf.xhk;
        end
        %st = [x1(i,j); y1(i,j)];
%         Logs{j}.sV_ekf(:,i)= st;
        % Compute squared error
        %Logs{j}.eV_ekf(:,i) = (jpdaf.config.TrackList{j}.TrackObj.pf.xhk(1:2,1) - st).*(jpdaf.config.TrackList{j}.TrackObj.pf.xhk(1:2,1) - st);
    end
    TrackNum_log(i) = TrackNum;
    if (ShowPlots)
        if(i==1 || rem(i,SkipFrames+1)==0)
            % Plot data
            cla(ax(1));
             % Flip the image upside down before showing it
            %imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            rectangle(ax(1),'Position',[8.75 4 0.75 3.5])
            rectangle(ax(1),'Position',[4 6 2 1])
            rectangle(ax(1),'Position',[1 3 2 2])
%             for j=1:TrueTracks
%                 h2 = plot(ax(1), x_true(1:i,j),y_true(1:i,j),'b.-','LineWidth',1);
%                 if j==2
%                     set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%                 end
%                 h2 = plot(ax(1), x_true(i,j),y_true(i,j),'bo','MarkerSize', 10);
%                 set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
%             end
            if(ShowData)
                data = plot(ax(1), DataList{i}(1,:),DataList{i}(2,:),'g*','MarkerSize', 5);
            end
            for j=1:TrackNum
                colour = 'r';
%                 if(j==2)
%                    colour = 'c';
%                 elseif (j==3)
%                    colour = 'm';
%                 end
                h4 = plot(ax(1), jpdaf.config.TrackList{j}.TrackObj.pf.xhk(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.xhk(2,:),strcat(colour,'.-'),'LineWidth',1);
                %h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
                c_mean = mean(jpdaf.config.TrackList{j}.TrackObj.pf.particles,2);
                c_cov = [std(jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.w')^2,0;0,std(jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),jpdaf.config.TrackList{j}.TrackObj.pf.w')^2];
                h2=plot_gaussian_ellipsoid(c_mean(1:2), c_cov, 1, [], ax(1));
                quiver(ax(1), c_mean(1),c_mean(2),30*c_mean(3)*cos(c_mean(4)),30*c_mean(3)*sin(c_mean(4)), 'Color', 'b', 'Linewidth', 1, 'ShowArrowHead','off', 'MaxHeadSize', 1)
                set(h2,'color',colour);
                set(h2,'LineWidth',1);
                %plot(ax(1),jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                text(ax(1),c_mean(1)+20,c_mean(2)-5,int2str(j));
                text(ax(1),c_mean(1)+20,c_mean(2)-15,num2str(jpdaf.config.TrackList{j}.TrackObj.pf.ExistProb, 2));
            end
                % set the y-axis back to normal.
            %set(ax(1),'ydir','normal');
            str = sprintf('Visualisation of tracking process');
            title(ax(1),str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
%            h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
%            set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
            axis(ax(1),V_bounds)
            %pause(0.01)
            if(Record)
                
                F(i) = getframe(ax(1));
            end
            
            % Plot PHD
            cla(ax(2), 'reset');
            [bandwidth,density,X,Y]=kde2d(myphd.config.particles(1:2,:)');
            %contour3(X,Y,density,50);
            surf(ax(2),X,Y,density);
            hold on;
            plot(ax(2), myphd.config.particles(1,:), myphd.config.particles(2,:), '.')
            hold on;
            plot(ax(2), myphd.config.z(1,:), myphd.config.z(2,:), 'y*');
            axis(ax(2),V_bounds);
            str = sprintf('Visualisation of PHD search track density');
            xlabel(ax(2),'X position (m)')
            ylabel(ax(2),'Y position (m)')
            zlabel(ax(2),'Intensity')
            title(ax(2),str)
            %pause(0.01)
            
            
        end
    end
end

F = F(2:end);
vidObj = VideoWriter(sprintf('phd_search.avi'));
vidObj.Quality = 100;
vidObj.FrameRate = 10;
open(vidObj);
writeVideo(vidObj, F);
close(vidObj);

% % Plot surveillance region
% figure('units','centimeters','position',[.1 .1 30 20])
% imagesc([min_x max_x], [min_y max_y], flipud(img));
% hold on
% rectangle('Position',[8.75 4 0.75 3], 'FaceColor', [0.8 0.8 0.8])
% text(8.8,6.8,'A_3', 'FontWeight', 'bold')
% rectangle('Position',[4 6 2 1], 'FaceColor', [0.8 0.8 0.8])
% text(4.05,6.8,'A_2', 'FontWeight', 'bold')
% rectangle('Position',[1 3 2 2], 'FaceColor', [0.8 0.8 0.8])
% text(1.05,4.8,'A_1', 'FontWeight', 'bold')
% h1 = plot(x1(:,1), y1(:,1), 'rx--', 'LineWidth', 2)
% h2 = plot(x_obstructed(:,1), y_obstructed(:,1), 'r--', 'LineWidth', 2)
% h3 = plot(x1(:,2), y1(:,2), 'bx--', 'LineWidth', 2)
% h4 = plot(x_obstructed(:,2), y_obstructed(:,2), 'b--', 'LineWidth', 2)
% h5 = plot(x1(:,3), y1(:,3), 'cx--', 'LineWidth', 2)
% h6 = plot(x_obstructed(:,3), y_obstructed(:,3), 'c--', 'LineWidth', 2)
% set(gca, 'ydir','normal');
% h_legend = legend([h1, h3, h5], {'Target 1', 'Target 2', 'Target3'}, 'Orientation','horizontal');
% %tit = "RMS Positional error (\sigma_{m}=";
% %tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
% %title(tit)
% ylabel("y-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
% xlabel("x-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
% set(gca,'FontSize',20)
% set(h_legend,'FontSize',16);
% axis([0 10 0 10])
% box on
% 
% % Plot measurement scan
% figure('units','centimeters','position',[.1 .1 30 20])
% imagesc([min_x max_x], [min_y max_y], flipud(img));
% hold on
% h1 = plot(DataList{1}(1,:),DataList{1}(2,:),'r*','MarkerSize', 10);
% h2 = plot(DataList{1}(1,1:3),DataList{1}(2,1:3),'k*','MarkerSize', 10);
% set(gca, 'ydir','normal');
% h_legend = legend([h2, h1], {'True measurements', 'Clutter measurements'}, 'Orientation','horizontal');
% %tit = "RMS Positional error (\sigma_{m}=";
% %tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
% %title(tit)
% ylabel("y-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
% xlabel("x-coordinate (m)", 'FontSize', 20, 'FontWeight', 'bold')
% set(gca,'FontSize',20)
% set(h_legend,'FontSize',16);
% axis([0 10 0 10])
% box on

% Plot number of track_estimates
figure('units','centimeters','position',[.1 .1 30 20])
hold on
TrueTrackNum_log = sum(~isnan(x1),2);
h1 = area(1:1184, TrueTrackNum_log)
h2 = plot(1:1184, permute(sim_res(2,2,:),[2 3 1]), 'r-','LineWidth', 3)
h3 = plot(1:1184, permute(sim_res(2,1,:),[2 3 1]), 'g-', 'LineWidth', 3)
plot(repmat(135,1,6), [0:1:5], 'y--')
text(135,0.1,'T_1\downarrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(175,1,6), [0:1:5], 'y--')
text(175,0.1,'T_1\uparrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(265,1,6), [0:1:5], 'y--')
text(265,0.1,'T_3\downarrow', 'FontWeight', 'bold', 'Color', 'yellow', 'HorizontalAlignment','right', 'FontSize', 12)
plot(repmat(306,1,6), [0:1:5], 'y--')
text(306,0.1,'T_3\uparrow', 'FontWeight', 'bold', 'Color', 'yellow', 'HorizontalAlignment','right', 'FontSize', 12)
plot(repmat(339,1,6), [0:1:5], 'y--')
text(339,0.1,'T_2\downarrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(380,1,6), [0:1:5], 'y--')
text(380,0.1,'T_2\uparrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(432,1,6), [0:1:5], 'y--')
text(432,0.1,'T_3\downarrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(471,1,6), [0:1:5], 'y--')
text(471,0.1,'T_3\uparrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(630,1,6), [0:1:5], 'y--')
text(630,0.1,'T_3\downarrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(729,1,6), [0:1:5], 'y--')
text(729,0.1,'T_2\downarrow', 'FontWeight', 'bold', 'Color', 'yellow', 'HorizontalAlignment','right', 'FontSize', 12)
plot(repmat(770,1,6), [0:1:5], 'y--')
text(770,0.1,'T_2\uparrow', 'FontWeight', 'bold', 'Color', 'yellow', 'HorizontalAlignment','right', 'FontSize', 12)
plot(repmat(788,1,6), [0:1:5], 'y--')
text(788,0.1,'T_3\uparrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(987,1,6), [0:1:5], 'y--')
text(987,0.1,'T_2\downarrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
plot(repmat(1046,1,6), [0:1:5], 'y--')
text(1046,0.1,'T_2\uparrow', 'FontWeight', 'bold', 'Color', 'yellow', 'FontSize', 12)
h_legend = legend([h1, h2, h3], {'True tracks', 'PDAF-LLR-TM', 'PHD-EP-TM'}, 'Orientation','horizontal');
tit = "Number of targets vs. Time (\lambda_{FA}V=50)";
%tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
title(tit)
ylabel("Number of tracks", 'FontSize', 20, 'FontWeight', 'bold')
xlabel("Time (s)", 'FontSize', 20, 'FontWeight', 'bold')
set(gca,'FontSize',20)
set(h_legend,'FontSize',16);
axis([0 1184 0 5])
box on

% Plot execution times
figure('units','centimeters','position',[.1 .1 30 20])
hold on
h1 = plot(lambda_list, exec_times(1,:), 'gx-');
h2 = plot(lambda_list, exec_times(2,:), 'rx-');
h_legend = legend([h2, h1], {'PDAF-LLR-TM', 'PHD-EP-TM'}, 'Orientation','horizontal');
tit = "Execution time vs. Clutter rate";
%tit = strcat(tit, sprintf("%0.2fm)",r_list(r_iter)));
title(tit)
ylabel("Time (s)", 'FontSize', 20, 'FontWeight', 'bold')
xlabel("\lambda_{FA}V", 'FontSize', 20, 'FontWeight', 'bold')
set(gca,'FontSize',20)
set(h_legend,'FontSize',16);
axis([0 100 0 700])
box on
%TrueTrackNum_log = zeros(size(x_true,1),1);
%for i = 1:size(x_true,2)
%end
