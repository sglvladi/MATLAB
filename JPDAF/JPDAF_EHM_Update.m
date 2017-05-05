%   JPDAF_EHM_Update.m                           Author: Lyudmil Vladimirov
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
%   Dependencies: buildEHMnet2.m 
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function TrackList = JPDAF_EHM_Update(TrackList,DataList,ValidationMatrix, bettaNTFA)

%% Initiate parameters
TrackNum    = size(TrackList,2); % Number of targets including FA
% alpha       = 0.3;      % log likelihood forget factor
PG          = 0.9;      % probability of Gating
PD          = 0.8;      % probability of Detection
GateLevel   = 5;
PointNum = size(ValidationMatrix,1);

%% Form clusters of tracks sharing measurements
clusters = {};
% Measurement Clustering
for i=1:TrackNum % Iterate over all tracks (including dummy)
    matched =[];
    temp_clust = find(ValidationMatrix(:,i))'; % Extract associated tracks
    
    % If measurement matched with any tracks
    if (~isempty(temp_clust))   
        % Check if matched tracks are members of any clusters
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
                clusters{1,match_ind}=union(cell2mat(clusters(1,match_ind)), cell2mat(clusters(1,match_ind+1)));
                clusters(:,match_ind+1)=[];
            end
            % Finally, join with associated track.
            clusters{1,match_ind}=union(cell2mat(clusters(1,match_ind)), temp_clust);
        else % If not matched with any cluster, then create a new one.
            clusters{1,size(clusters,2)+1} = temp_clust;
        end
    end
end

% Dummy 1 cluster generation
% for i=1:TrackNum % Iterate over all tracks (including dummy)
%     matched =[];
%     temp_clust = find(ValidationMatrix(:,i))'; % Extract associated tracks
%     clusters{1,1}=union(cell2mat(clusters(1,1)), temp_clust);
%     clusters{1,1} = clusters{1,1}(:)';
% end

% Build ClusterList
ClusterList = [];
ClusterObj.MeasIndList = [];
ClusterObj.TrackIndList = [];
for c=1:size(clusters,2)
    ClusterList{c} = ClusterObj;
    ClusterList{c}.MeasIndList = clusters{1,c}(:)';
    for i = 1:size(ClusterList{c}.MeasIndList,2)
        ClusterList{c}.TrackIndList = union(ClusterList{c}.TrackIndList, find(ValidationMatrix(ClusterList{c}.MeasIndList(i)',:)));
    end
end

%% Compute Association Likelihoods 
Li = zeros(PointNum, TrackNum);
%Li(:,1) = ones(size(Li,1), 1)*bettaNTFA*(1-PD*PG); 
for i=1:PointNum
    for j=1:TrackNum
        z = DataList(:,i);
        z_pred = TrackList{j}.TrackObj.z_pred;
        S = TrackList{j}.TrackObj.S;
        Li(i,j) = mvnpdf(z, z_pred, S)*PD;
    end
end
%Li(PointNum+1,:) = ones(1, size(Li,2))*bettaNTFA*(1-PD*PG);


%% Create Hypothesis net for each cluster
%NodeList = [];
NetList = [];
NetObj.NodeList = [];
NetObj.EdgeList = []; 

% Hypothesis net for each cluster
for c=1:size(ClusterList,2)
    Cluster = ClusterList{c};
    ClustMeasIndList = Cluster.MeasIndList;
    ClustTrackIndList = Cluster.TrackIndList;
    NetList{c} = buildEHMnet2(ValidationMatrix(ClustMeasIndList', ClustTrackIndList), [ones(size(Li(ClustMeasIndList', ClustTrackIndList),1), 1)*bettaNTFA*(1-PD*PG),Li(ClustMeasIndList', ClustTrackIndList)]);
end

%NetObj = buildEHMnet(ValidationMatrix, Li);%(ClustMeasIndList', ClustTrackIndList)

%% Compute weights and update each track
for i=1:TrackNum,
    
    % Get the index of the cluster which track belongs to
    for j=1:size(ClusterList,2)
        Cluster = ClusterList{j};
        if (ismember(i, Cluster.TrackIndList)~=0)
            cluster_id = j;
            break;
        end
    end
    
    % Get the EHM Net relating to that cluster
    NetObj = NetList{cluster_id};
    %NetObj = NetList{1};
    
   
    % Extract prediction information (Computed from ObservationAssociation.m) 
    x_pred       = TrackList{i}.TrackObj.x_pred;    % predicted mean
    P_pred       = TrackList{i}.TrackObj.P_pred;    % predicted covariance
    z_pred       = TrackList{i}.TrackObj.z_pred;    % predicted measurement
    S            = TrackList{i}.TrackObj.S;         % predicted measurement covariance
    Sinv         = inv(S);
    DataInd      = find(ValidationMatrix(:,i))';    % Associated measurements

    % extract measurements
    z = DataList(:,DataInd);
    [ObsDim,ValidDataPointNum] = size(z);
    
    % Kalman Filter
    innov_err   = z - z_pred(:,ones(1,ValidDataPointNum)); % error (innovation) for each sample
    W           = TrackList{i}.TrackObj.P12*Sinv;                    % Kalman gain matrix
    %Pc          = (eye(size(F,1)) - W*H)*P_pred;
    Pc          = P_pred - W*S*W';

    
    %------------------------------------------------
    % compute association probabilities (non-parametric JPDA)

    % separate memory for likelihood ratios and association prababilities
    %Li_k    = zeros(1, ValidDataPointNum+1);   
    betta   = zeros(1,ValidDataPointNum); 
    
    % Compute likelihood ratios
    ClustMeasIndList=[];
    for j=1:size(DataInd,2)
        ClustMeasIndList(j) = find(ClusterList{cluster_id}.MeasIndList==DataInd(j));
    end    
    ClustTrackInd = find(ClusterList{cluster_id}.TrackIndList==i)+1; % T1 is the false alarm
    
    % Extract betta for target
    betta = NetObj.betta(ClustMeasIndList,ClustTrackInd)'

    %------------------------------------------------
    % update
    tot_innov_err    = innov_err*betta(1:ValidDataPointNum)';
    Pgag    = W*((innov_err.*betta(ones(ObsDim,1),1:ValidDataPointNum))*innov_err' - tot_innov_err*tot_innov_err')*W';
    Pnew    = (1-sum(betta,2))*P_pred + (sum(betta,2))*Pc + Pgag;
    xnew    = x_pred + W*tot_innov_err;

    
   % Update track
    TrackList{i}.TrackObj.x     = xnew;
    TrackList{i}.TrackObj.P     = Pnew;
    
end;    % track loop
