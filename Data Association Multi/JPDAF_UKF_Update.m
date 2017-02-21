function [TrackList, betta_overall] = JPDAF_UKF_Update(TrackList,DataList,ValidationMatrix, bettaNTFA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Track_Update - performs track update
% Input:
%   TrackList    - List of all target tracks at time k(TrackObj's)
%   DataList     - List of all measurements at time k
%   EKF          - Kalman Filter Structure
% Output:
%   TrackList    - Updated TrackList
%
% [1] Y. Bar-Shalom, F. Daum, and J. Huang, 
%     "Probabilistic Data Association Filter", 
%     IEEE Control Systems Magazine, December 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
TrackNum    = size(TrackList,2);
% alpha       = 0.3;      % log likelihood forget factor
PG          = 0.9;      % probability of Gating
PD          = 0.8;      % probability of Detection
GateLevel   = 5;
PointNum = size(ValidationMatrix,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form clusters of tracks sharing measurements
clusters = {}

UnclusteredTracks = [1:TrackNum];
for i=1:PointNum % Iterate over all measurements
    matched =[];
    temp_clust = find(ValidationMatrix(i,:))'; % Extract associated tracks
    UnclusteredTracks = setdiff(UnclusteredTracks,temp_clust);
    if (~isempty(temp_clust))   % If measurement matched with any tracks
        % Check if matched tracks are members of any clusters
        for j=1:size(clusters,2)
            a = ismember(temp_clust, cell2mat(clusters(1,j)));
            if (ismember(1,a)~=0)
                matched = [matched, j]; % Store matched cluster ids
            end   
        end
        if(size(matched,2)==1) % If only matched with a single cluster, join.
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

if(~isempty(UnclusteredTracks))
    for i=1:size(UnclusteredTracks)
        clusters{end+1} = UnclusteredTracks(i);
    end
end

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
%clusters(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop throught all tracks
betta_overall = zeros(TrackNum, PointNum);
for i=1:TrackNum,
    for j=1:size(clusters,2)
        cell2mat(clusters(1,j))
        a = ismember(i, cell2mat(clusters(1,j)));
        if (a==1)
            cluster_id = j;
        end   
    end
    %------------------------------------------------
    % extract prediction
    x_pred       = TrackList{i}.TrackObj.x_pred;
    P_pred       = TrackList{i}.TrackObj.P_pred;
    z_pred       = TrackList{i}.TrackObj.z_pred;
    S            = TrackList{i}.TrackObj.S;
    Sinv         = inv(S);
    %F            = TrackList{i}.TrackObj.F;
    %H            = TrackList{i}.TrackObj.H;
    DataInd      = find(ValidationMatrix(:,i))';

    %------------------------------------------------
    % extract measurements
    z = DataList(:,DataInd);
    [ObsDim,ValidDataPointNum] = size(z);
    
    %----------------------------------
    % basic kalman filter
    innov_err   = z - z_pred(:,ones(1,ValidDataPointNum)); % error (innovation) for each sample
    W           = TrackList{i}.TrackObj.P12*Sinv;                    % Kalman gain matrix
    %Pc          = (eye(size(F,1)) - W*H)*P_pred;
    Pc          = P_pred - W*S*W';

    
    %------------------------------------------------
    % compute association probabilities (as described in [1] for a non-parametric PDA)
    
    %C   = pi; % volume of the 2-dimensional unit hypersphere     
    % V_k = C*GateLevel^(ObsDim/2)*det(S)^(1/2)   % volume of the validation region 
    
    % separate memory for likelihood ratios and association prababilities
    Li_k    = zeros(1, ValidDataPointNum+1);   
    betta   = zeros(1,ValidDataPointNum+1); 
    i
    DataInd
    Li
    % Compute likelihood ratios
    for j=1:size(DataInd,2)
        betta(j) = getAssocProb(DataInd(j),i, NetList{cluster_id}.HypMap, NetList{cluster_id}.HypProb);
        betta_overall(i,DataInd(j)) = betta(j);
        %Li_k(j) =  mvnpdf(z(:,j), z_pred, S)*PD/(ValidDataPointNum/V_k); % Likelihood ratio of measurement i
    end
    
    
    
    % Compute association probabilities
    betta(ValidDataPointNum+1) = getAssocProb(PointNum+1,i, NetList{cluster_id}.HypMap, NetList{cluster_id}.HypProb);
    betta
    % Back-up code
    % ========================================================================>
    % loglik  = sum((innov_err'*Sinv).*innov_err',2); % volume computation for the vector
    % % The next is OK
    % betta(1:ValidDataPointNum) = exp(-.5*loglik); %a
    % betta(ValidDataPointNum+1) = (1-PG*PD)/PD*2*ValidDataPointNum/GateLevel*sqrt(det(S)); % by Tracking and Data Association
    % betta(vlen+1) = vlen*(1-PG*PD)/PD/(pi*gamma*sqrt(det(S))); % by Multitarget-Multisensor Tracking
    % betta   = betta./sum(betta);
    % ========================================================================/
    
    %------------------------------------------------
    % update
    tot_innov_err    = innov_err*betta(1:ValidDataPointNum)';
    Pgag    = W*((innov_err.*betta(ones(ObsDim,1),1:ValidDataPointNum))*innov_err' - tot_innov_err*tot_innov_err')*W';
    Pnew    = betta(ValidDataPointNum+1)*P_pred + (1-betta(ValidDataPointNum+1))*Pc + Pgag;
        
    xnew    = x_pred + W*tot_innov_err;

    
   % This is not a real PDAF update
    TrackList{i}.TrackObj.x     = xnew;
    TrackList{i}.TrackObj.P     = Pnew;
    
end;    % track loop
