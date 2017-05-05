%   JPDAF_UKF_Update.m                           Author: Lyudmil Vladimirov
%   ======================================================================>
%   Functionality: 
%       Compute association weights (betta) and perform track
%       update for each target.
%   
%   Input: 
%       TrackList    - List of all target tracks at time k(TrackObj's)
%       DataList     - List of all measurements at time k
%       ValidationMatrix    - Matrix containing all possible measurement to
%                             track associations.
%                             (Output from ObservationAssociation.m)
%       bettaNTFA    - New track/False alarm density (assumed to be same)
%       
%   
%   Dependencies: buildHyp.m, getAssocProb.m 
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [TrackList, betta_overall] = JPDAF_UKF_Update(TrackList,DataList,ValidationMatrix, bettaNTFA)

%% Initiate parameters
TrackNum    = size(TrackList,2);    % Number of targets
% alpha       = 0.3;      % log likelihood forget factor
PG          = 0.9;      % probability of Gating
PD          = 0.8;      % probability of Detection
GateLevel   = 5;        % Gate Level 
PointNum = size(ValidationMatrix,1);    % Number of measurements


%% Form clusters of tracks sharing measurements
clusters = {};

UnclusteredTracks = [1:TrackNum];
for i=1:PointNum % Iterate over all measurements
    matched =[];
    temp_clust = find(ValidationMatrix(i,:))'; % Extract associated tracks
    UnclusteredTracks = setdiff(UnclusteredTracks,temp_clust);
    
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

% If there still exist Unclustered tracks
if(~isempty(UnclusteredTracks))
    for i=1:size(UnclusteredTracks)
        clusters{end+1} = UnclusteredTracks(i); % Create new cluster for each 
    end
end

%% Compute Likelihood table
Li = zeros(PointNum+1, TrackNum);
for i=1:PointNum
    for j=1:TrackNum
        z = DataList(:,i);
        z_pred = TrackList{j}.TrackObj.z_pred;
        S = TrackList{j}.TrackObj.S;
        Li(i,j) = mvnpdf(z, z_pred, S)*PD;
    end
end
Li(PointNum+1,:) = ones(1, size(Li,2))*bettaNTFA*(1-PD*PG);


%% Create Hypothesis tree for each cluster
%NodeList = [];
NetList = [];
HypTree.HypMap = [];
HypTree.HypProb = [];
% Hypothesis net for each cluster
for c=1:size(clusters,2)
    NetList{c} = HypTree;
    ClustTracks = clusters{1,c};
    NetList{c}.HypMap = containers.Map('KeyType','double', 'ValueType','any');
    NetList{c}.HypMap(1) = [];
    NetList{c}.HypProb = [bettaNTFA];
    [NetList{c}.HypMap, NetList{c}.HypProb] = buildHyp(NetList{c}.HypMap, NetList{c}.HypProb, Li, ClustTracks, ValidationMatrix, [], bettaNTFA);
end

%% Compute weights and update each track
betta_overall = zeros(TrackNum, PointNum);
for i=1:TrackNum
    
    % Identify the cluster where track i belongs to
    for j=1:size(clusters,2)
        cell2mat(clusters(1,j));
        a = ismember(i, cell2mat(clusters(1,j)));
        if (a==1)
            cluster_id = j;
        end   
    end
    
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
    Li_k    = zeros(1, ValidDataPointNum+1);   
    betta   = zeros(1,ValidDataPointNum+1); 

    % Compute likelihood ratios
    for j=1:size(DataInd,2)
        betta(j) = getAssocProb(DataInd(j),i, NetList{cluster_id}.HypMap, NetList{cluster_id}.HypProb);
        betta_overall(i,DataInd(j)) = betta(j);
        %Li_k(j) =  mvnpdf(z(:,j), z_pred, S)*PD/(ValidDataPointNum/V_k); % Likelihood ratio of measurement i
    end
    
    % Compute betta for FA measurement
    betta(ValidDataPointNum+1) = getAssocProb(PointNum+1,i, NetList{cluster_id}.HypMap, NetList{cluster_id}.HypProb);
    
    %------------------------------------------------
    % update
    tot_innov_err    = innov_err*betta(1:ValidDataPointNum)';
    Pgag    = W*((innov_err.*betta(ones(ObsDim,1),1:ValidDataPointNum))*innov_err' - tot_innov_err*tot_innov_err')*W';
    Pnew    = betta(ValidDataPointNum+1)*P_pred + (1-betta(ValidDataPointNum+1))*Pc + Pgag;
        
    xnew    = x_pred + W*tot_innov_err;

    
   % Update track
    TrackList{i}.TrackObj.x     = xnew;
    TrackList{i}.TrackObj.P     = Pnew;
    
end;    % track loop
