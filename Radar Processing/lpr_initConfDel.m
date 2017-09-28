%% Plot settings
ShowPlots = 1;
SkipFrames = 0;

%% Initiate PF parameters
nx = 4;      % number of state dims
nu = 4;      % size of the vector of process noise
nv = 2;      % size of the vector of observation noise
q  = 0.5;   % process noise density (std)
r  = 2;    % observation noise density (std)
lambdaV = 5; % mean number of clutter points 
% Prior PDF generator
gen_x0_cch = @(Np) mvnrnd(repmat([0,0,0,0],Np,1),diag([q^2, q^2, 100, 100]));
% Process equation x[k] = sys(k, x[k-1], u[k]);
sys_cch = @(k, xkm1, uk) [xkm1(1,:)+2*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+2*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ 2*uk(:,3)'; xkm1(4,:) + 2*uk(:,4)'];
% PDF of process noise generator function
gen_sys_noise_cch = @(u) mvnrnd(zeros(size(u,2), nu), diag([0,0,q^2,0.3^2])); 
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
pf.Np              = 5000;                 % number of particles
pf.particles       = zeros(5, pf.Np); % particles
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

%% Set TrackNum
TrackNum = 0;
TrueTracks = 3;

%% Generate DataList                       (TrackNum, x_true, y_true, R, R_clutter, lambdaV, Iter)                   
%[DataList,x1,y1] = gen_obs_cluttered_multi2(TrueTracks, x_true, y_true, r, 2, lambdaV, 1);

%% Get GroundTruth
for i=1:TrueTracks
    GroundTruth{i} = [x_true(:,i), y_true(:,i)]; % ith target's GroundTruth
end

%% Initiate TrackList
% for i=1:TrackNum,
%     pf.gen_x0 = @(Np) mvnrnd(repmat([GroundTruth{i}(1,1),GroundTruth{i}(1,2),0,0],Np,1),diag([q^2, q^2, 1, 1]));
%     pf.ExistProb = 0.8;
%     TrackList{i}.TrackObj = ParticleFilterMin2(pf);
% end;

%% Initiate JPDAF parameters
Par = [];
Par.Filter = ParticleFilterMin2(pf);
Par.DataList = DataList{1}(:,:);
Par.GroundTruth = GroundTruth;
Par.TrackList = [];
Par.PD = 0.8;
Par.PG = 0.998;
Par.GateLevel = 5;
Par.Pbirth = 0.001;
Par.Pdeath = 0.4;
Par.SimIter = 1000;

%% Assign PHD parameter values
par.k               = 1;                                                    % initial iteration number
par.Np              = 10000;                                                % number of particles
par.resampling_strategy = 'systematic_resampling';                          % resampling strategy
par.birth_strategy = 'mixture';                                           %  
par.sys = @(k, xkm1, uk) [xkm1(1,:)+2*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+2*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)']; % CH model
par.gen_x0 = @(Np) [10*rand(Np,1),10*rand(Np,1), mvnrnd(zeros(Np,1), 2*sigma_v^2), 2*pi*rand(Np,1)]; % Uniform position and heading, Gaussian speed
par.particles = par.gen_x0(par.Np)';                                        % Generate inital particles as per gen_x0
par.w = repmat(1/par.Np, par.Np, 1)';                                       % Uniform weights
par.likelihood = @(k, yk, xk) mvnpdf(yk, xk, cov_v);                        % Likelihood model p(y|x)
par.obs_model = @(xk) [xk(1,:); xk(2,:)];                                   % Observation model (no noise)
par.sys_noise = gen_sys_noise_cch;                                          % System noise
par.Pbirth = 0.1;                                                           % Birth probability
par.Pdeath = 0.005;                                                         % Death probability
par.J_k = 5000;                                                             % Number of birth particles (births by "expansion")
par.PD = 0.9;                                                               % Probability of detection
par.lambda = lambdaV/(10^2);                                                % Mean clutter per unit area
par.Np_conf = pf.Np;                                                        % Number of particles for confirmed tracks
par.P_conf = 0.9;                                                           % Confirmation probability
par.type = 'search';                                                        % Search PHD filter

myphd = SMC_PHD(par);

%% Instantiate JPDAF
jpdaf = JPDAF(Par);
Par.pdaf = 1;
jpdaf_init = JPDAF(Par);

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
    
    % 1) Predict the confirmed tracks
    try
        jpdaf.Predict();
    catch
        error('This');
    end
    % 2) Update the confirmed track
    jpdaf.Update();
    
    % 3) Perform Track Management (Score-based)
    [jpdaf, jpdaf_init] = Track_InitConfDel(jpdaf, jpdaf_init);
    
    TrackNum = size(jpdaf.config.TrackList,2);
    
    exec_time = exec_time + toc;
    
    %store Logs
    for j=1:TrackNum
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
%             for j=1:TrackNum
%                 h2 = plot(Logs{j}.sV_ekf(1,1:i),Logs{j}.sV_ekf(2,1:i),'b.-','LineWidth',1);
%                 if j==2
%                     set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%                 end
%                 h2 = plot(Logs{j}.sV_ekf(1,i),Logs{j}.sV_ekf(2,i),'bo','MarkerSize', 10);
%                 set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
%             end
            h2 = plot(ax(1), DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
            for j=1:TrackNum
                colour = 'b';
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
                %plot(jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                text(c_mean(1)+20,c_mean(2)-10,int2str(j));
            end
            
            for j=1:jpdaf_init.config.TrackNum
                colour = 'r';
                h4 = plot(ax(1), jpdaf_init.config.TrackList{j}.TrackObj.pf.xhk(1,:),jpdaf_init.config.TrackList{j}.TrackObj.pf.xhk(2,:),strcat(colour,'.-'),'LineWidth',1);
                %h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
                c_mean = mean(jpdaf_init.config.TrackList{j}.TrackObj.pf.particles,2);
                c_cov = [std(jpdaf_init.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf_init.config.TrackList{j}.TrackObj.pf.w')^2,0;0,std(jpdaf_init.config.TrackList{j}.TrackObj.pf.particles(2,:),jpdaf_init.config.TrackList{j}.TrackObj.pf.w')^2];
                h2=plot_gaussian_ellipsoid(c_mean(1:2), c_cov, 1, [], ax(1));
                set(h2,'color',colour);
                set(h2,'LineWidth',1);
                %plot(jpdaf.config.TrackList{j}.TrackObj.pf.particles(1,:),jpdaf.config.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
                % set the y-axis back to normal.
            set(ax(1),'ydir','normal');
            str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
            title(str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
%            h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
%            set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
            axis(ax(1),[-2000 200 -3000 3000])
            %axis(ax(1),[-1800 -1600 -200 0])
%             % Plot PHD
%             cla(ax(2), 'reset');
%             [bandwidth,density,X,Y]=kde2d(myphd.config.particles(1:2,:)');
%             %contour3(X,Y,density,50);
%             surf(ax(2),X,Y,density);
%             hold on;
%             plot(ax(2), myphd.config.particles(1,:), myphd.config.particles(2,:), '.')
%             hold on;
%             plot(ax(2), myphd.config.z(1,:), myphd.config.z(2,:), 'y*');
%             axis(ax(2), [0 10 0 10 0 100]);
             pause(0.1)
        end
    end
end

figure
plot(1:1184, TrackNum_log)
hold on
%TrueTrackNum_log = zeros(size(x_true,1),1);
%for i = 1:size(x_true,2)
TrueTrackNum_log = sum(~isnan(x1),2);
%end
plot(1:1184, TrueTrackNum_log)