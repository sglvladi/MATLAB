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

function TrackList = JPDAF_TO_EHM_PF_Update(TrackList,DataList, ValidationMatrix, bettaNTFA, search)

%% Initiate parameters
TrackNum    = size(TrackList,2); % Number of targets including FA
% alpha       = 0.3;      % log likelihood forget factor
PG          = 0.9933;      % probability of Gating
PD          = 0.7;      % probability of Detection
GateLevel   = 10;
PointNum = size(ValidationMatrix,1);
clustering  = 1;

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
        %else if no measurements exists for the given track    
%         else
%             % Check if an empty measurement cluster already exists 
%             for j=1:size(clusters,2)
%                 if(isempty(cell2mat(clusters(1,j))))
%                     matched = [matched, j]; % Store matched cluster ids
%                 end   
%             end
%             % If the above is true.
%             if(size(matched,2)==1) 
%                 %Nothing needs to be done
%                 %clusters{1,matched(1)}=union(cell2mat(clusters(1,matched(1))), temp_clust); 
%             else
%                 clusters{1,size(clusters,2)+1} = temp_clust; %create new cluster
%             end
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
    ClusterList{c}.MeasIndList = unique(clusters{1,c}(:)');
    
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
%Li(:,1) = ones(size(Li,1), 1)*bettaNTFA*(1-PD*PG); 
%parfor i=1:PointNum
    for j=1:TrackNum
        % Get valid measurement data indices
        ValidDataInd = find(TrackList{j}.TrackObj.pf.Validation_matrix(1,:));

        z = DataList(:,:);
        try
            z_pred = TrackList{j}.TrackObj.pf.z_pred;
        catch
            disp('error');
        end
        %S = TrackList{j}.TrackObj.pf.S;
        
        %Li_i(:,j) = mvnpdf(z', z_pred', S)*PD;
        %[Li(:,j),TrackList{j}.TrackObj.pf.Li']
        for i=1:size(ValidDataInd,2)
            Li(ValidDataInd(1,i),j) = TrackList{j}.TrackObj.pf.Li(1,i)'*PD*PG;
        end
        %[Li_i(:,j),Li(:,j)]
        
        %a=5;
    end
%end
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
    try
        ClustLi = [ones(size(Li(ClustMeasIndList', ClustTrackIndList)',1), 1)*bettaNTFA*(1-PD*PG),Li(ClustMeasIndList', ClustTrackIndList)']; 
    catch
        disp('error');
    end
    % Augment ValidationMatrix
%ValidationMatrix = [ones(1, PointNum); ValidationMatrix']';
    NetList{c} = buildEHMnet_trans([ones(size(ValidationMatrix(ClustMeasIndList', ClustTrackIndList)',1),1), ValidationMatrix(ClustMeasIndList', ClustTrackIndList)'], ClustLi);
    %dummy = buildEHMnet2([ones(1, size(ValidationMatrix(ClustMeasIndList', ClustTrackIndList),2));ValidationMatrix(ClustMeasIndList', ClustTrackIndList)'], ClustLi');
end

%NetObj = buildEHMnet(ValidationMatrix, Li);%(ClustMeasIndList', ClustTrackIndList)

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



        % Extract prediction information (Computed from ObservationAssociation.m) 
        %x_pred       = TrackList{i}.TrackObj.x_pred;    % predicted mean
        %P_pred       = TrackList{i}.TrackObj.P_pred;    % predicted covariance
        %z_pred       = TrackList{i}.TrackObj.z_pred;    % predicted measurement
        %S            = TrackList{i}.TrackObj.S;         % predicted measurement covariance
        %Sinv         = inv(S);
        DataInd      = find(ValidationMatrix(:,i))';    % Associated measurements

        % extract measurements
        z = DataList(:,DataInd);
        %[ObsDim,ValidDataPointNum] = size(z);

        % Compute likelihood ratios
        ClustMeasIndList=[];
        for j=1:size(DataInd,2)
           ClustMeasIndList(j) = unique(find(ClusterList{cluster_id}.MeasIndList==DataInd(j)));
        end    
        ClustTrackInd = find(ClusterList{cluster_id}.TrackIndList==i); % T1 is the false alarm

        % Extract betta for target
%         if(isempty(NetObj.betta_trans(ClustTrackInd, find(NetObj.betta_trans(ClustTrackInd, :)))))
%             disp('error');
%         end
        TrackList{i}.TrackObj.pf.ValidDataInd;
        TrackList{i}.TrackObj.pf.betta = [NetObj.betta(ClustTrackInd,1), NetObj.betta(ClustTrackInd, ClustMeasIndList+1)];
    else
        % Else if target was not matched with any clusters, it means it was
        % also not matched with any measurements and thus only the "dummy"
        % measurement association is possible (i.e. betta = [1]);
        TrackList{i}.TrackObj.pf.betta = 1;
    end
     
    %------------------------------------------------
    % update
    if(search)
        TrackList{i}.TrackObj.pf = TrackList{i}.TrackObj.UpdateSearch(TrackList{i}.TrackObj.pf);
    else
        TrackList{i}.TrackObj.pf = TrackList{i}.TrackObj.UpdateMulti(TrackList{i}.TrackObj.pf);
    end
    
end;    % track loop
