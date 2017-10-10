classdef JPDAFX <handle
% =====================================================================================
% Parameters:
% Par: structure with the following fields
%       
%       * Variables
%       -------------------
%       .DataList         = all available observations 
%                           (optional if GroundTruth is provided)
%       .TrackList        = A list of initiated targets (optional)
%       .TrackNum         = Number of targets to be generated
%                           (only applies when DataList is not provided)
%       .lambda           = False alarm density
%       .PD               = Probability of detection
%       .PG               = PG
%       .GateLevel        = GateLevel (Based on mahal. distance)
%       .InitDelMethod    = Track management method
%       .Pbirth           = Probability of new target birth
%       .Pdeath           = Probability of target death
%       .SimIter          = Number of timesteps to allow simulation for

    properties
        config
    end
    
    methods
        function obj = JPDAFX(config)
            
            % TrackList ~~~~~~>
            if isfield(config,'TrackList')
                config.TrackNum = size(config.TrackList,2);
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            % Validate .PD, .GateLevel, .Pbirth, .Pdeath, .SimIter ~~~~~~>
            if ~isfield(config,'PD') || ~isfield(config,'PG') || ~isfield(config,'GateLevel')
                error('One of the following has not been provide: PD, PG, GateLevel!');
            end
            % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>
            
            obj.config = config;
      
        end
        
        function Predict(obj)
            
            % Get number of available tracks and observations
            obj.config.TrackNum = size(obj.config.TrackList,2);
            obj.config.ObsNum   = size(obj.config.Y,2);
            
            % Validation matix and volume
            ValidationMatrix = zeros(obj.config.TrackNum, obj.config.ObsNum); % (Nt x Nm)
            LikelihoodMatrix = zeros(obj.config.TrackNum, obj.config.ObsNum); % (Nt x Nm)
            V_k = 0;
            
            if(~isempty(obj.config.TrackList))
                
                % Predict and construct Validation and Likelihood matrices
                for t = 1:obj.config.TrackNum
                    
                    TrackObj = obj.config.TrackList{t}.TrackObj;
                    TrackObj.config.k = obj.config.k;
                    
                    % Predict
                    TrackObj.Predict();
                    
                    % Extract predicted measurement and innovation covariance
                    if(isa(TrackObj,'ParticleFilterX'))
                        trans_parts = TrackObj.obs_model.obs(TrackObj.config.k, TrackObj.config.particles, TrackObj.obs_model.obs_noise(TrackObj.config.k,TrackObj.config.Np));
                        y_pred      = sum(repmat(TrackObj.config.w,size(trans_parts,1),1).*trans_parts,2);
                        S           = weightedcov(trans_parts',TrackObj.config.w');
                    elseif(isa(TrackObj,'KalmanFilterX'))
                        y_pred  = TrackObj.config.y_pred;
                        S       = TrackObj.config.S;
                    else
                        error('Class not defined!');
                    end
                    
                    % Perform Gating
                    switch numel(y_pred)
                        case 1 
                            C = 2;
                        case 2
                            C = pi;
                        case 3
                            C = 4*pi/3;
                        otherwise
                            error('[JPDAF] Gating has only been implemented for observations of up to 3 dimensions!');
                    end
                    V_k = V_k + C*obj.config.GateLevel^(numel(y_pred)/2)*det(S)^(1/2);    
                    ValidationMatrix(t,:) = obj.mahalDist(obj.config.Y, y_pred, S, 2) < obj.config.GateLevel;
                    
                    % Extract valid measurements
                    ValidDataInd = find(ValidationMatrix(t,:));
                    ValidY = obj.config.Y(:,ValidDataInd);
                    
                    % Update Likelihood matrix
                    if(isa(TrackObj,'ParticleFilterX'))
                        TrackObj.config.LikelihoodMatrix = TrackObj.obs_model.eval_likelihood(TrackObj.config.k, ValidY, TrackObj.config.particles);
                        LikelihoodMatrix(t, ValidDataInd) = obj.config.PD*obj.config.PG*sum(TrackObj.config.LikelihoodMatrix,2)'/TrackObj.config.Np;
                    elseif(isa(TrackObj,'KalmanFilterX'))
                        LikelihoodMatrix(t, ValidDataInd) = obj.config.PD*obj.config.PG*TrackObj.obs_model.eval_likelihood(TrackObj.config.k, ValidY, TrackObj.config.x_pred)';
                    end
                    TrackObj.config.y = ValidY;
                    obj.config.TrackList{t}.TrackObj = TrackObj;
                end
                
                obj.config.ValidationMatrix = ValidationMatrix;
                obj.config.LikelihoodMatrix = LikelihoodMatrix;
                
                % Compute New Track/False Alarm density
                obj.config.bettaNTFA = sum(obj.config.ValidationMatrix(:))/V_k;
                if(obj.config.bettaNTFA==0)
                    obj.config.bettaNTFA=1;
                end
                
                % Get all clusters
                UnassocTracks = obj.FormClusters();
                
                % Allocate memory for betta and fill in weights for unassociated tracks
                betta = zeros(obj.config.TrackNum, obj.config.ObsNum+1);
                for t = 1:numel(UnassocTracks)
                    betta(t,1) = 1;
                end

                % Create Hypothesis net for each cluster and populate betta
                NetList = cell(1,size(obj.config.ClusterList,2));
                for c=1:size(obj.config.ClusterList,2)
                    Cluster = obj.config.ClusterList{c};
                    ClustMeasIndList = Cluster.MeasIndList;
                    ClustTrackIndList = Cluster.TrackIndList;
                    ClustLi = zeros(numel(ClustTrackIndList), numel(ClustMeasIndList)+1);
                    ClustLi(:,1) = ones(numel(ClustTrackIndList),1)*obj.config.bettaNTFA*(1-obj.config.PD*obj.config.PG);
                    ClustLi(:,2:end) = obj.config.LikelihoodMatrix(ClustTrackIndList,ClustMeasIndList);
                    ClustVm = zeros(numel(ClustTrackIndList), numel(ClustMeasIndList)+1);
                    ClustVm(:,1) = ones(numel(ClustTrackIndList),1);
                    ClustVm(:,2:end) = obj.config.ValidationMatrix(ClustTrackIndList, ClustMeasIndList);
                    NetList{c} = buildEHMnet_trans(ClustVm, ClustLi);
                    betta(ClustTrackIndList, [1, ClustMeasIndList+1]) = NetList{c}.betta;
                end
                
                obj.config.NetList = NetList;
                obj.config.betta = betta;
            else
                fprintf('No tracks where found. Skipping JPDAF Predict step...\n');
                obj.config.ValidationMatrix = zeros(1, size(obj.config.DataList,2));
                obj.config.bettaNTFA = 0;
                obj.config.betta = -1; % Set betta to -1
            end
        end
        
        function Update(obj)
            if(~isempty(obj.config.TrackList))
                % Initiate parameters
                TrackNum    = size(obj.config.TrackList,2); % Number of targets including FA

                % Compute weights and update each track
                for t=1:TrackNum
                    ValidDataInd    = find(obj.config.ValidationMatrix(t,:));    % Associated measurements
                    % Only update if there exist available measurements
                    %if ~isempty(ValidDataInd)
                    betta = [obj.config.betta(t,1), obj.config.betta(t,ValidDataInd+1)];% [NetObj.betta(ClustTrackInd,1), NetObj.betta(ClustTrackInd, ClustMeasIndList+1)];
                    obj.config.TrackList{t}.TrackObj.UpdateMulti(betta);
                    %end
                end    
            else
                fprintf('No tracks where found. Skipping JPDAF Update step...\n');
            end
        end
        
        function UnassocTracks = FormClusters(obj)
            % Initiate parameters
            TrackNum    = size(obj.config.TrackList,2); % Number of targets including FA
            ValidationMatrix = obj.config.ValidationMatrix;
            clustering  = 1;

            % Form clusters of tracks sharing measurements
            clusters = {};
            UnassocTracks = [];
            if(clustering)
                if(isfield(obj.config, 'pdaf'))
                    % Do nothing
                else
                    % Measurement Clustering
                    for i=1:TrackNum % Iterate over all tracks 
                        matched =[];
                        temp_clust = find(ValidationMatrix(i,:)); % Extract associated measurements

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
                        else
                            UnassocTracks = [UnassocTracks i];
                        end
                    end
                end
            else
                % Measurement Clustering
                for i=1:TrackNum % Iterate over all tracks

                    temp_clust = find(ValidationMatrix(i,:)); % Extract associated tracks

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
            if(isfield(obj.config, 'pdaf'))
                for i=1:TrackNum
                    ClusterList{i} = ClusterObj;
                    ClusterList{i}.MeasIndList = find(ValidationMatrix(i,:));
                    ClusterList{i}.TrackIndList = i;
                end
            else
                for c=1:size(clusters,2)
                    ClusterList{c} = ClusterObj;
                    ClusterList{c}.MeasIndList = unique(clusters{1,c}(:)');

                    % If we are currently processing the cluster of unassociated tracks 
                    if(isempty(ClusterList{c}.MeasIndList))
                        ClusterList{c}.TrackIndList = unique(union(ClusterList{c}.TrackIndList, find(all(ValidationMatrix==0))));
                    else
                        for i = 1:size(ClusterList{c}.MeasIndList,2) 
                            ClusterList{c}.TrackIndList = unique(union(ClusterList{c}.TrackIndList, find(ValidationMatrix(:,ClusterList{c}.MeasIndList(i)))));
                        end
                    end
                end
            end
            obj.config.ClusterList = ClusterList;
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