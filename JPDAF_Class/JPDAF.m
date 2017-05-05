classdef JPDAF
% =====================================================================================
% Parameters:
% Par: structure with the following fields
%       
%       * Variables
%       -------------------
%       .GroundTruth      = ground truth data (optional)
%       .DataGenerator    = Data generator class instance
%                           (optional if DataList is provided)
%       .DataList         = all available observations 
%                           (optional if GroundTruth is provided)
%       .TrackList        = A list of initiated targets (optional)
%       .TrackNum         = Number of targets to be generated
%                           (only applies when DataList is not provided)
%       .Filter           = Underlying filter class instance to be used
%                           (KF, EKF, UKF or PF - only used if TrackList is not provided)
%       .lambda           = False alarm density
%       .PD               = Probability of detection
%       .PG               = PG
%       .GateLevel        = GateLevel (Based on mahal. distance)
%       .InitDelMethod    = Track management method
%       .Pbirth           = Probability of new target birth
%       .Pdeath           = Probability of target death
%       .SimIter          = Number of timesteps to allow simulation for

    properties
        Par
    end
    
    methods
        function obj = JPDAF(prop)
            % Validate .Filter ~~~~~~~~~~~~~~~~~~~~~~>
            if ~isfield(prop,'Filter')&&~isfield(prop,'TrackList')
                error('Base Filter class instance (Par.Filter) has not been provided.. Please instantiate the desired filter (KF, EKF, UKF or PF) and include it as an argument! \n');             
            elseif ~isfield(prop,'TrackList')
                if(isa(prop.Filter,'ParticleFilterMin2'))
                    prop.FilterType = 'PF';
                elseif(isa(prop.Filter,'KalmanFilter_new'))
                    prop.FilterType = 'KF';
                elseif(isa(prop.Filter,'EKalmanFilter'))
                    prop.FilterType = 'UKF';
                elseif(isa(prop.Filter,'UKalmanFilter'))
                    prop.FilterType = 'EKF';
                else
                    error('Base Filter class instance (Par.Filter) is invalid.. Please instantiate the desired filter (KF, EKF, UKF or PF) and include it as an argument! \n');
                end
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            % Validate .DataList, .GroundTruth ~~~~~~>
            if ~isfield(prop,'DataList') && ~isfield(prop,'GroundTruth')
                error('No DataList and no Ground Truth have been supplied. please provide one of the two in order to proceed.\nNOTE: if only Ground Truth is supplied, then a DataGenerator instance needs to also be provided.');
            elseif ~isfield(prop,'DataList') && ~isfield(prop,'DataGenerator')
                error('If only Ground Truth is supplied, then a DataGenerator instance needs to also be provided such that a simulated DataList can be produced.');    
            elseif ~isfield(prop,'DataList')
                prop.DataList = prop.DataGenerator.genDataList();
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            % TrackList ~~~~~~>
            if isfield(prop,'TrackList')
                prop.TrackNum = size(prop.TrackList,2);
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            
            % Validate .PD, .GateLevel, .Pbirth, .Pdeath, .SimIter ~~~~~~>
            if ~isfield(prop,'PD') || ~isfield(prop,'PG') || ~isfield(prop,'GateLevel') || ~isfield(prop,'Pbirth') || ~isfield(prop,'Pdeath') || ~isfield(prop,'SimIter')
                error('One of the following has not been provide: PD, PG, GateLevel, Pbirth, Pdeath, SimIter!');
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            obj.Par = prop;
      
        end
        
        function Par = Predict(~, Par)
            if(isa(Par.TrackList{1}.TrackObj,'ParticleFilterMin2'))
                Par.DataList( :, ~any(Par.DataList,1) ) = []; 
                %[TrackList, ValidationMatrix, bettaNTFA] = Observation_Association(TrackList, tempDataList, ekf);
                %TrackList = JPDAF_EKF_Update(TrackList, DataList(:,:,i), ValidationMatrix', bettaNTFA);
                %[TrackList, betta] = JPDAF_UKF_Update(TrackList, tempDataList, ValidationMatrix', bettaNTFA);
                %TrackList = Track_InitConfDel(TrackList,tempDataList,ValidationMatrix',bettaNTFA, betta);
                Par.ValidationMatrix = zeros(Par.TrackNum, size(Par.DataList,2)); 
                tot_gate_area = 0;
                for t = 1:Par.TrackNum
                    %Par.TrackList{t}.TrackObj.pf.k = i;
                    Par.TrackList{t}.TrackObj.pf.z = Par.DataList;
                    Par.TrackList{t}.TrackObj.pf = Par.TrackList{t}.TrackObj.PredictMulti(Par.TrackList{t}.TrackObj.pf);
                    Par.ValidationMatrix(t,:) = Par.TrackList{t}.TrackObj.pf.Validation_matrix;
                    tot_gate_area = tot_gate_area + Par.TrackList{t}.TrackObj.pf.V_k;
                end
                % Compute New Track/False Alarm density
                Par.bettaNTFA = sum(Par.ValidationMatrix(:))/tot_gate_area;
            end
        end
        
        function Par = Update(obj, Par)
            if(isa(Par.TrackList{1}.TrackObj,'ParticleFilterMin2'))
                [Par.TrackList] = obj.PF_Update(Par);
            end
        end
        
        %   PF_Update                           Author: Lyudmil Vladimirov
        %   ======================================================================>
        %   Functionality: 
        %       Compute association weights (betta) and perform track
        %       update for each target, using EHM.
        %   
        %   Input: 
        %       TrackList    - List of all target tracks at time k(TrackObj's)
        %       DataList     - List of all measurements at time k
        %       ValidationMatrix    - Matrix containing all possible measurement to
        %                             track associations.
        %                             (Output from ObservationAssociation.m)
        %       bettaNTFA    - New track/False alarm density (assumed to be same)
        %   
        %   Output:
        %       TrackList    - Updated list of all target tracks
        %   
        %   Dependencies: buildEHMnet_Fast.m 
        %   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        function TrackList = PF_Update(~,Par)
            %% Initiate parameters
            TrackNum    = size(Par.TrackList,2); % Number of targets including FA
            % alpha       = 0.3;      % log likelihood forget factor
            PG          = Par.PG;      % probability of Gating
            PD          = Par.PD;      % probability of Detection
            GateLevel   = Par.GateLevel;
            bettaNTFA    = Par.bettaNTFA;
            ValidationMatrix = Par.ValidationMatrix';
            PointNum = size(Par.ValidationMatrix,1);
            TrackList = Par.TrackList;
            clustering  = 1;

            if(~isempty(ValidationMatrix)&&(~isempty(find(sum(ValidationMatrix,1)==0,1))))
                disp('Check');
            end
            %% Form clusters of tracks sharing measurements
            clusters = {};
            if(clustering)
                % Measurement Clustering
                for i=1:TrackNum % Iterate over all tracks 
                    matched =[];
                    temp_clust = find(ValidationMatrix(:,i))'; % Extract associated measurements

                    % If track matched with any measurements
                    if (~isempty(temp_clust))   
                        % Check if matched measurements are members of any clusters
                        for j=1:size(clusters,2)
                            a = ismember(temp_clust, cell2mat(clusters(1,j)));
                            if (ismember(1,a)~=0)
                                matched = [matched, j]; % Store matched cluster ids
                            end   
                        end

                        % If only matched with a single cluster, join.
                        if(size(matched,2)==1) 
                            clusters{1,matched(1)}=union(cell2mat(clusters(1,matched(1))), temp_clust);
                        elseif (size(matched,2)>1) % If matched with more that one clusters
                            matched = sort(matched); % Sort cluster ids
                            % Start from last cluster, joining each one with the previous
                            %   and removing the former.  
                            for match_ind = size(matched,2)-1:-1:1
                                clusters{1,matched(match_ind)}=union(cell2mat(clusters(1,matched(match_ind))), cell2mat(clusters(1,matched(match_ind+1))));
                                clusters(:,matched(match_ind+1))=[];
                            end
                            % Finally, join with associated track.
                            clusters{1,matched(match_ind)}=union(cell2mat(clusters(1,matched(match_ind))), temp_clust);
                        else % If not matched with any cluster, then create a new one.
                            clusters{1,size(clusters,2)+1} = temp_clust;
                        end
                     end
                end
            else
                % Measurement Clustering
                for i=1:TrackNum % Iterate over all tracks (including dummy)

                    temp_clust = find(ValidationMatrix(:,i))'; % Extract associated tracks

                    % If measurement matched with any tracks
                    if (~isempty(temp_clust))
                        if(~isempty(clusters))
                            clusters{1,1}= union(clusters{1,1}, temp_clust);
                        else
                           clusters{1,1}= temp_clust; 
                        end
                    end
                end
            end

            % Build ClusterList
            ClusterList = [];
            ClusterObj.MeasIndList = [];
            ClusterObj.TrackIndList = [];
            for c=1:size(clusters,2)
                ClusterList{c} = ClusterObj;
                ClusterList{c}.MeasIndList = fast_unique(clusters{1,c}(:)');

                % If we are currently processing the cluster of unassociated tracks 
                if(isempty(ClusterList{c}.MeasIndList))
                    ClusterList{c}.TrackIndList = unique(union(ClusterList{c}.TrackIndList, find(all(ValidationMatrix==0))));
                else
                    for i = 1:size(ClusterList{c}.MeasIndList,2) 
                        ClusterList{c}.TrackIndList = unique(union(ClusterList{c}.TrackIndList, find(ValidationMatrix(ClusterList{c}.MeasIndList(i)',:))));
                    end
                end
            end

            %% Compute Association Likelihoods 
            Li = zeros(PointNum, TrackNum);
            for j=1:TrackNum
                % Get valid measurement data indices
                ValidDataInd = find(TrackList{j}.TrackObj.pf.Validation_matrix(1,:));
                z = Par.DataList(:,:);
                try
                    z_pred = TrackList{j}.TrackObj.pf.z_pred;
                catch
                    disp('error');
                end
                for i=1:size(ValidDataInd,2)
                    Li(ValidDataInd(1,i),j) = TrackList{j}.TrackObj.pf.Li(1,i)'*PD;
                end
            end


            %% Create Hypothesis net for each cluster
            NetList = [];
            NetObj.NodeList = [];
            NetObj.EdgeList = []; 

            % Hypothesis net for each cluster
            for c=1:size(ClusterList,2)
                Cluster = ClusterList{c};
                ClustMeasIndList = Cluster.MeasIndList;
                ClustTrackIndList = Cluster.TrackIndList;
                try
                    ClustLi = [ones(size(Li(ClustMeasIndList', ClustTrackIndList),1), 1)*bettaNTFA*(1-PD*PG),Li(ClustMeasIndList', ClustTrackIndList)]; 
                catch
                    disp('error');
                end
                NetList{c} = buildEHMnet_fast(ValidationMatrix(ClustMeasIndList', ClustTrackIndList), ClustLi);
            end

            %% Compute weights and update each track
            for i=1:TrackNum,

                cluster_id = 0;
                % Get the index of the cluster which track belongs to
                for j=1:size(ClusterList,2)
                    Cluster = ClusterList{j};
                    if (ismember(i, Cluster.TrackIndList)~=0)
                        cluster_id = j;
                        break;
                    end
                end

                % If target has been matched with a cluster
                %  then extract it's association prob. matrix
                if(cluster_id~=0)
                    try
                        % Get the EHM Net relating to that cluster
                        NetObj = NetList{cluster_id};
                        %NetObj = NetList{1};
                    catch
                        disp('this');
                    end

                    DataInd      = find(ValidationMatrix(:,i))';    % Associated measurements

                    % extract measurements
                    z = Par.DataList(:,DataInd);

                    % Compute likelihood ratios
                    ClustMeasIndList=[];
                    for j=1:size(DataInd,2)
                        ClustMeasIndList(j) = fast_unique(find(ClusterList{cluster_id}.MeasIndList==DataInd(j)));
                    end    
                    ClustTrackInd = find(ClusterList{cluster_id}.TrackIndList==i); % T1 is the false alarm

                    % Extract betta for target
                    if(isempty(NetObj.betta_trans(ClustTrackInd, find(NetObj.betta_trans(ClustTrackInd, :)))))
                        disp('error');
                    end
                    TrackList{i}.TrackObj.pf.betta = NetObj.betta_trans(ClustTrackInd, find(NetObj.betta_trans(ClustTrackInd, :)));
                else
                    % Else if target was not matched with any clusters, it means it was
                    % also not matched with any measurements and thus only the "dummy"
                    % measurement association is possible (i.e. betta = [1]);
                    TrackList{i}.TrackObj.pf.betta = 1;
                end

                %------------------------------------------------
                % update
                TrackList{i}.TrackObj.pf = TrackList{i}.TrackObj.UpdateMulti(TrackList{i}.TrackObj.pf);
            end    % track loop
        end
        
       function Par =  TrackInitConfDel(~,Par)
            for t = 1:Par.TrackNum
                if(Par.TrackList{t}.TrackObj.pf.ExistProb>0.1)
                    TrackList{end} = Par.TrackList{t};
                end
            end
            TrackNum = size(TrackList,2);

            invalidDataInd = find((sum(Par.ValidationMatrix,1)==0));
            % Process search track
            if(Par.pf_search.pf.ExistProb>0.9)
                disp('Search Track Exist prob:');
                Par.pf_search.pf.z = Par.DataList(:, invalidDataInd);
                Par.pf_search.pf = Par.pf_search.PredictSearch(Par.pf_search.pf);
                Par.pf_search.pf = Par.pf_search.UpdateSearch(Par.pf_search.pf);

                if(Par.pf_search.pf.ExistProb>0.9)
                    % Promote new track
                    TrackNum = TrackNum + 1;
                    TrackList{TrackNum}.TrackObj = Par.pf_search;

                    % Create new PF search track
                    nx = 4;      % number of state dims
                    nu = 4;      % size of the vector of process noise
                    nv = 2;      % size of the vector of observation noise
                    q  = 0.01;   % process noise density (std)
                    r  = 0.3;    % observation noise density (std)
                    % Process equation x[k] = sys(k, x[k-1], u[k]);
                    sys_cch = @(k, xkm1, uk) [xkm1(1,:)+1*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+1*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:)+ uk(:,3)'; xkm1(4,:) + uk(:,4)'];
                    % PDF of process noise generator function
                    gen_sys_noise_cch = @(u) mvnrnd(zeros(size(u,2), nu), diag([0,0,q^2,0.16^2])); 
                    % Observation equation y[k] = obs(k, x[k], v[k]);
                    obs = @(k, xk, vk) [xk(1)+vk(1); xk(2)+vk(2)];                  % (returns column vector)
                    % PDF of observation noise and noise generator function
                    sigma_v = r;
                    cov_v = sigma_v^2*eye(nv);
                    p_obs_noise   = @(v) mvnpdf(v, zeros(1, nv), cov_v);
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
                    pf.obs = p_yk_given_xk;
                    pf.obs_model = @(xk) [xk(1,:); xk(2,:)];
                    pf.R = cov_v;
                    pf.clutter_flag = 1;
                    pf.multi_flag = 0;
                    pf.sys_noise = gen_sys_noise_cch;
                    pf.gen_x0 = @(Np) [10*rand(Np,1),10*rand(Np,1), mvnrnd(zeros(Np,1), 2*sigma_v^2), 2*pi*rand(Np,1)];
                    %pf.xhk = [s.x_init(1,i),s.x_init(2,i),0,0]';
                    pf.ExistProb = 0.5;
                    Par.pf_search = ParticleFilterMin2(pf);

                    disp('Promoted one track');
                end
            else
                disp('Search Track Exist prob:');
                Par.pf_search.pf.z = Par.DataList(:, invalidDataInd);
                Par.pf_search.pf = Par.pf_search.PredictSearch(Par.pf_search.pf);
                Par.pf_search.pf = Par.pf_search.UpdateSearch(Par.pf_search.pf);;
                if(Par.pf_search.pf.ExistProb<0.1)
                    % Reset the search track
                    pf.gen_x0 = @(Np) [10*rand(Np,1),10*rand(Np,1), mvnrnd(zeros(Np,1), 2*sigma_v^2), 2*pi*rand(Np,1)];
                    %pf.xhk = [s.x_init(1,i),s.x_init(2,i),0,0]';
                    pf.ExistProb = 0.5;
                    pf_search = ParticleFilterMin2(pf);
                    pf_search.pf.multi_flag = 0;
                end

            end
        end
    end
    
    methods (Static)
        function D=mahalDist(x, m, C, use_log)
        % p=gaussian_prob(x, m, C, use_log)
        %
        % Evaluate the multi-variate density with mean vector m and covariance
        % matrix C for the input vector x.
        % Vectorized version: Here X is a matrix of column vectors, and p is 
        % a vector of probabilities for each vector.

            if nargin<4, use_log = 0; end

            d   = length(m);

            if size(x,1)~=d
               x=x';
            end
            N       = size(x,2);

            m       = m(:);
            M       = m*ones(1,N);
            denom   = (2*pi)^(d/2)*sqrt(abs(det(C)));
            invC    = inv(C);
            mahal   = sum(((x-M)'*invC).*(x-M)',2);   % Chris Bregler's trick

            switch use_log,
            case 2,
              D     = mahal;
            case 1,
              D     = -0.5*mahal - log(denom);
            case 0,
              numer = exp(-0.5*mahal);
              D     = numer/denom;
            otherwise
                error('Unsupported log type')
            end
        end
    end
end