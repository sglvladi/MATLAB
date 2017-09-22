function [jpdaf, jpdaf_init] = Track_InitConfDel_MN(jpdaf, jpdaf_init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track_Maintenance - Performs track initiation, confirmation and deletion 
% Input:
%   TrackList        - List of Tracks
%   ValidationMatrix - Matrix showing all valid data associations
% Output:
%   TrackList    - Updated list of Tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Confirmed tracks
    
    TrackList = jpdaf.config.TrackList;
    DataList = jpdaf.config.DataList;
    if (isfield(jpdaf.config, 'ValidationMatrix'))
        ValidationMatrix = jpdaf.config.ValidationMatrix;
        bettaNTFA = jpdaf.config.bettaNTFA;
        betta = jpdaf.config.betta;
    else
        ValidationMatrix = zeros(size(TrackList,2), size(DataList,2));
        bettaNTFA = 0;
        betta = [];
    end

    TrackNum    = size(TrackList,2);
    PointNum    = size(DataList,2);
    ObsDim      = size(DataList,1);
    Vmax = 0.4;     
    
    % High and low thresholds for confirmation and deletion
    %  of tentative tracks
    ConfThresh  = 4;       
    DelThresh   = 2;
    
    % if the ValidationMatrix is empty (signals overall initiation step)
    if numel(ValidationMatrix)==0
        UnassocMeasInd = [1:PointNum];
    else
        % Get indices of all rows (measurements) where there exist no possible association
        UnassocMeasInd = find(all(ValidationMatrix==0,1));
    end

    % Update LPR and status of existing existing tracks
    %ConfirmedTracks = find(cell2mat(cellfun(@(sas)sas.Status, TrackList.TrackObj, 'uni', false))=='Confirmed');
    DeletedTracks = [];
    for i=1:TrackNum
        TrackInd = i;
        PossibleAssocMeas = find(ValidationMatrix(TrackInd,:));

        % Update LPR
        if isempty(PossibleAssocMeas)
            TrackList{TrackInd}.TrackObj.pf.Det = 0;
            TrackList{TrackInd}.TrackObj.pf.Miss = TrackList{TrackInd}.TrackObj.pf.Miss + 1;
            %fprintf('LPR Track %d = %f\n', TrackInd, TrackList{TrackInd}.TrackObj.pf.LPR);
        else
            TrackList{TrackInd}.TrackObj.pf.Det = TrackList{TrackInd}.TrackObj.pf.Det+1;
            TrackList{TrackInd}.TrackObj.pf.Miss = 0;
        end

        % Check against thresholds
        if (TrackList{TrackInd}.TrackObj.pf.Miss > DelThresh)
            DeletedTracks(end+1) = TrackInd;
            %TrackList{TrackInd} = []; % Delete Track
        end  
    end
    
    DeletedTracks = sort(DeletedTracks);
    for i=1:numel(DeletedTracks)
        TrackList(DeletedTracks(numel(DeletedTracks)+1-i)) = [];
    end
%     if (~isempty(TrackList))
%         % Remove cells of deleted tracks from cell array
%         TrackList(~cellfun('isempty',TrackList));
%     end
    
    jpdaf.config.TrackList = TrackList;
    jpdaf.config.TrackNum = size(TrackList,2);
    
    
    %% Tentative tracks 
    try
        jpdaf_init.config.DataList = DataList(:,UnassocMeasInd); % Only use unassociated measurements
    catch
        error('safg');
    end
    % 1) Predict the confirmed tracks
    jpdaf_init.Predict();
    % 2) Update the confirmed track
    jpdaf_init.Update();
    
    TrackList = jpdaf_init.config.TrackList;
    DataList = jpdaf_init.config.DataList;
    if (isfield(jpdaf_init.config, 'ValidationMatrix'))
        ValidationMatrix = jpdaf_init.config.ValidationMatrix;
        bettaNTFA = jpdaf_init.config.bettaNTFA;
        betta = jpdaf_init.config.betta;
    else
        ValidationMatrix = zeros(size(TrackList,2), size(DataList,2));
        bettaNTFA = 0;
        betta = [];
    end

    TrackNum    = size(TrackList,2);
    PointNum    = size(DataList,2);
    ObsDim      = size(DataList,1);
    Vmax = 0.4;
    
    % Kalman Parameters
    n=4;      %number of state
    q=0.01;    %std of process 
    r=0.1;    %std of measurement
    s.Q=[1^3/3, 0, 1^2/2, 0;  0, 1^3/3, 0, 1^2/2; 1^2/2, 0, 1, 0; 0, 1^2/2, 0, 1]*10*q^2; % covariance of process
    s.R=r^2*eye(n/2);        % covariance of measurement  
    s.sys=(@(x)[x(1)+ x(3); x(2)+x(4); x(3); x(4)]);  % assuming measurements arrive 1 per sec
    s.obs=@(x)[x(1);x(2)];                               % measurement equation                                % initial state
    s.x_init = [];
    s.P_init = [];
    
    % Initiate PF parameters
    nx = 4;      % number of state dims
    nu = 4;      % size of the vector of process noise
    nv = 2;      % size of the vector of observation noise
    q  = 0.01;   % process noise density (std)
    r  = 0.1;    % observation noise density (std)
    
    lambdaV = 50; % mean number of clutter points 
    % Prior PDF generator
    gen_x0_cch = @(Np) mvnrnd(repmat([0,0,0,0],Np,1),diag([q^2, q^2, 100, 100]));
    % Process equation x[k] = sys(k, x[k-1], u[k]);
    sys_cch = @(k, xkm1, uk) [xkm1(1,:)+1*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+1*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)'];
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
    pf.ExistProb = 0.5;
    
    % Two-point difference initiation
    %for i = 1:TrackNum
    %    s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); DataList(1,i,2)-DataList(1,i,1); DataList(2,i,2)-DataList(2,i,1)]; %initial state
    %end
    %s.P_init = [q^2, 0, q^2, 0;
    %            0, q^2, 0, q^2;
    %            q^2, 0, 2*q^2, 0;
    %            0, q^2, 0, 2*q^2];                               % initial state covraiance

    % Single-point initiatio
    
    % if the ValidationMatrix is empty (signals overall initiation step)
    if numel(ValidationMatrix)==0
        UnassocMeasInd = [1:PointNum];
    else
        % Get indices of all rows (measurements) where there exist no possible association
        UnassocMeasInd = find(all(ValidationMatrix==0,1));
    end

    % Update LPR and status of existing existing tracks
    %ConfirmedTracks = find(cell2mat(cellfun(@(sas)sas.Status, TrackList.TrackObj, 'uni', false))=='Confirmed');
    DeletedTracks = [];
    for i=1:TrackNum
        TrackInd = i;
        PossibleAssocMeas = find(ValidationMatrix(TrackInd,:));

        if isempty(PossibleAssocMeas)
            TrackList{TrackInd}.TrackObj.pf.Det = 0;
            TrackList{TrackInd}.TrackObj.pf.Miss = TrackList{TrackInd}.TrackObj.pf.Miss + 1;
            %fprintf('LPR Track %d = %f\n', TrackInd, TrackList{TrackInd}.TrackObj.pf.LPR);
        else
            TrackList{TrackInd}.TrackObj.pf.Det = TrackList{TrackInd}.TrackObj.pf.Det+1;
            TrackList{TrackInd}.TrackObj.pf.Miss = 0;
        end

        % Check against thresholds
        if (TrackList{TrackInd}.TrackObj.pf.Status==0) % Tentative tracks
            if TrackList{TrackInd}.TrackObj.pf.Det>ConfThresh
                TrackList{TrackInd}.TrackObj.pf.Status = 1; % Confirm track
                jpdaf.config.TrackList{end+1} = TrackList{TrackInd};
                DeletedTracks(end+1) = TrackInd;
            elseif (TrackList{TrackInd}.TrackObj.pf.Miss>DelThresh)
                DeletedTracks(end+1) = TrackInd;
                %TrackList{TrackInd} = []; % Delete Track
            end
        end
    end
    
    DeletedTracks = sort(DeletedTracks);
    for i=1:numel(DeletedTracks)
        try
            TrackList(DeletedTracks(numel(DeletedTracks)+1-i)) = [];
        catch
            error('Here');
        end
    end
%     if (~isempty(TrackList))
%         % Remove cells of deleted tracks from cell array
%         TrackList(~cellfun('isempty',TrackList));
%     end
        
    % Initiate new tracks
    for i=1:numel(UnassocMeasInd)
        %warning('Initiating track');
        MeasInd = UnassocMeasInd(i);
        pf.gen_x0 = @(Np) mvnrnd(repmat([DataList(1,MeasInd), DataList(2,MeasInd), 0, 0],Np,1),diag([q^2, q^2, 0.3^2, 4*pi^2]));
        %s.x = [DataList(1,MeasInd); DataList(2,MeasInd); 0; 0]; %initial state
        %s.P = diag([q^2, q^2, (Vmax^2/3), (Vmax^2/3)]);
        pf.Status = 0;
        pf.Det = 1;
        pf.Miss = 0;
        %pf.LPR = 1;
        %pf.LPR_max = pf.LPR;
        TrackList{end+1}.TrackObj = ParticleFilterMin2(pf); 
    end
    jpdaf_init.config.TrackList = TrackList;
    jpdaf_init.config.TrackNum = size(TrackList,2);
    jpdaf.config.TrackNum = size(jpdaf.config.TrackList,2);
end