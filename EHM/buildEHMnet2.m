%   buildEHMnet2.m                              Author: Lyudmil Vladimirov
%   ======================================================================>
%   Functionality: 
%       Build EHM net and compute respective association probabilities (betta).
%       (normally executed for each cluster)
%
%   Input:
%       ValidationMatrix    - Matrix (m x n) containing all possible measurement
%                             to track associations (for cluster of interest).
%                             (Output from ObservationAssociation.m)
%       Li                  - Matrix (m x n) containing association likelihoods
%       (m: number of measurements in cluster)
%       (n: number of associated tracks, including dummy track)

%   Output:
%       Structure NetObj:
%           #NetObj.NodeList - List of all Nodes (NodeObj) contained in net
%           #NetObj.EdgeList - Cell Matrix (n x n, n being the total number of
%                              nodes), where cell (i,j) contains the tracks
%                              contained in edge from parent node i, to
%                              child node j.
%           #NetObj.ValidationMatrix - Internal Validation matrix for
%                                      cluster.
%           #NetObj.betta    - Matrix (n x m) containing Association  
%                              probabilites computed for all m tracks and
%                              n measurements.
%           #Betta for the transposed problem (m x n) can be computed as:
%               betta_trans = [ones(1,size(betta,2)-1)-sum(betta(:,2:end),1); betta(:,2:end)] 
%  
%   IMPORTANT REMINDER TO SELF:
%     When computing the remainders for a node, we always look at the remaining
%     association events in the !NEXT! layer and then get the difference
%     between these and any entries in the node's MeasIndList.
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function NetObj = buildEHMnet2(ValidationMatrix, Li)

% Get number of measurements/layers
PointNum = size(ValidationMatrix,1); 
LayerNum = PointNum; % Layer 1 is root layer

% Augment ValidationMatrix
%ValidationMatrix = [ones(1, PointNum); ValidationMatrix']';

% Get number of Tracks
TrackNum = size(ValidationMatrix,2);

% Define Node Object
NodeObj.MeasInd = 0; % Measurement index
NodeObj.TrackIndList = []; % List of associated tracks
NodeObj.ParentIndList = []; % List of Parent nodes
NodeObj.ChildIndList = [];  % List of Child nodes
NodeObj.Remainders = [ 1:TrackNum ]'; % Remaining tracks

% Define Net Object
NetObj.NodeList = [];
NetObj.EdgeList = [];
NetObj.ValidationMatrix = ValidationMatrix;
NetObj.Li = Li;

% Create Root Node
NetObj.NodeList{1} = NodeObj;

% Build the net
% ==============>
% For every measurement/layer
for j=1:LayerNum

    % Get index list (L_jm1_Ind) of nodes in previous layer (j-1)
    L_jm1_Ind = find(cell2mat(cellfun(@(sas)sas.MeasInd, NetObj.NodeList, 'uni', false))==j-1);

    % Get indeces of all associated targets for measurement j
    T_j = [find(ValidationMatrix(j,:))]; % i=0 for false alarm

    % For every node in L_jm1
    for i_jm1 = 1:size(L_jm1_Ind,2)

        % Index of parent node
        ParentInd = L_jm1_Ind(i_jm1);

        % Get all targets to consider
        T_jm1 = union(1,setdiff(T_j,NetObj.NodeList{ParentInd}.TrackIndList)); 

        % For every track in T_jm1
        for i=1:size(T_jm1,2)

            % Get the track index
            TrackInd = T_jm1(i);

            % Init Match flag
            match_flag = 0;

            dfg = cellfun(@(sas)sas.MeasInd, NetObj.NodeList, 'uni', false );
            dfg = cell2mat(dfg);
            % Get index list L_j of nodes in current layer (j)
            L_j_Ind = find(dfg==j);

            % if current Layer is empty
            if isempty(L_j_Ind)
                % Index of child node
                ChildInd = size(NetObj.NodeList,2)+1;

                % Create child node
                NetObj.NodeList{ChildInd} = NodeObj;
                NetObj.NodeList{ChildInd}.MeasInd = j;
                NetObj.NodeList{ChildInd}.TrackIndList = union(NetObj.NodeList{L_jm1_Ind(i_jm1)}.TrackIndList(:,:), T_jm1(i));
                NetObj.NodeList{ChildInd}.TrackIndList = NetObj.NodeList{ChildInd}.TrackIndList(:)';
                NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);

                % Create edge from parent to child
                NetObj.EdgeList{ParentInd, ChildInd} =[TrackInd];

                % Compute remainders
                T_rem_j = []; % T_j+1:mk[N^(j)_(i_j)]
                for j_sub = j+1:LayerNum
                    T_rem_j = union(T_rem_j, find(ValidationMatrix(j_sub,:))');
                end
                NetObj.NodeList{ChildInd}.Remainders = setdiff(T_rem_j,setdiff(NetObj.NodeList{ChildInd}.TrackIndList,1)); 
            else
         
                % Compute remainders (j-1)
                T_rem_jm1_ti = []; %T_j+1:mk[N^(j-1)_(i_j-1),t_i]
                for j_sub = j+1:LayerNum
                    T_rem_jm1_ti = union(T_rem_jm1_ti, find(ValidationMatrix(j_sub,:))');
                end
                R_jm1 = setdiff(T_rem_jm1_ti,setdiff(union(NetObj.NodeList{ParentInd}.TrackIndList, TrackInd),1));
                R_jm1 = R_jm1(:); % Enforce that R_jm1 is a column vector

                % For all nodes in current layer
                for i_j=1:size(L_j_Ind,2)
                    ChildInd = L_j_Ind(i_j);

                    % Compute remainders (j)
                    %T_rem_j = []; % T_j+1:mk[N^(j)_(i_j)]
                    %for j_sub = j+1:LayerNum
                    %    T_rem_j = find(ValidationMatrix(j+1:LayerNum,:))';
                    %end
                    R_j = NetObj.NodeList{ChildInd}.Remainders;

                    % If the node's list of remainders (R_j) is equal to R_jm1
                    if (isequal(R_jm1,R_j)) 
                        % Simply add a new edge and update parent/child
                        %   relationships
                        NetObj.NodeList{ChildInd}.TrackIndList = union(NetObj.NodeList{ChildInd}.TrackIndList,TrackInd);
                        if size(NetObj.EdgeList,1)<ParentInd
                            NetObj.EdgeList{ParentInd, ChildInd} = TrackInd;
                        else
                            NetObj.EdgeList{ParentInd, ChildInd} = union(NetObj.EdgeList{ParentInd, ChildInd}, TrackInd);
                        end
                        NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                        NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);
                        NetObj.NodeList{ChildInd}.TrackIndList = union(NetObj.NodeList{ChildInd}.TrackIndList, NetObj.NodeList{ParentInd}.TrackIndList);
                        
                        % set match_flag
                        match_flag = 1;
                        break;
                    end
                end

                if match_flag==0
                    % Index of child node
                    ChildInd = size(NetObj.NodeList,2)+1;

                    % Create child node
                    NetObj.NodeList{ChildInd} = NodeObj;
                    NetObj.NodeList{ChildInd}.MeasInd = j;
                    NetObj.NodeList{ChildInd}.TrackIndList = union(NetObj.NodeList{L_jm1_Ind(i_jm1)}.TrackIndList(:,:), T_jm1(i));
                    NetObj.NodeList{ChildInd}.TrackIndList = NetObj.NodeList{ChildInd}.TrackIndList(:)';
                    NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                    NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);

                    % Create edge from parent to child
                    NetObj.EdgeList{ParentInd, ChildInd} =[TrackInd];

                    % Compute remainders
                    T_rem_j = []; % T_j+1:mk[N^(j)_(i_j)]
                    for j_sub = j+1:LayerNum
                        T_rem_j = union(T_rem_j, find(ValidationMatrix(j_sub,:))');
                    end
                    R_j = setdiff(T_rem_j,setdiff(NetObj.NodeList{ChildInd}.TrackIndList,1));
                    NetObj.NodeList{ChildInd}.Remainders = R_j; 
                end
            end

        end
    end
end 

% Calculate the vector P_D
p_D = zeros(size(NetObj.NodeList,2),1);
p_D(1,1) = 1;
% for every node
for NodeInd=2:size(NetObj.NodeList,2)
    p_D_temp = 0;
    Node = NetObj.NodeList{NodeInd};
    
    % for every parent of the node
    for i = 1:size(Node.ParentIndList,2)
        
        ParentInd = Node.ParentIndList(i);
        p_D_m1 = p_D(ParentInd,1);
        
        % for every track in the parent-child edge
        TrackEdgeList = cell2mat(NetObj.EdgeList(ParentInd,NodeInd));
        for t = 1:size(TrackEdgeList,2)
            TrackInd = TrackEdgeList(t);
            p_D_temp = p_D_temp + Li(Node.MeasInd, TrackInd)*p_D_m1;
        end
    end
    p_D(NodeInd,1) = p_D_temp; 
end
NetObj.p_D = p_D; 
% Calculate the vector P_U
p_U = zeros(size(NetObj.NodeList,2),1);
p_U(end,1) = 1;
% for every node starting from the leaf
for i=1:size(NetObj.NodeList,2)-1
    NodeInd = (size(NetObj.NodeList,2)) - i;
    p_U_temp = 0;
    Node = NetObj.NodeList{NodeInd};
    
    % for every child of the node
    for j = 1:size(Node.ChildIndList,2)
        ChildInd = Node.ChildIndList(j);
        ChildNode = NetObj.NodeList{ChildInd};
        p_U_p1 = p_U(ChildInd,1);
        
        % for every track in the parent-child edge
        TrackEdgeList = cell2mat(NetObj.EdgeList(NodeInd, ChildInd));
        for t = 1:size(TrackEdgeList,2)
            TrackInd = TrackEdgeList(t);
            p_U_temp = p_U_temp + Li(ChildNode.MeasInd, TrackInd)*p_U_p1;
        end
    end
    p_U(NodeInd,1) = p_U_temp; 
end
NetObj.p_U = p_U;
% Compute P_DT matrix
p_DT = zeros(TrackNum, size(NetObj.NodeList,2));
for NodeInd=1:size(NetObj.NodeList,2)
    Node = NetObj.NodeList{NodeInd};
    p_DT_temp = 1;
     
    for TrackInd = 1:TrackNum
        
        % Valid parents, where TrackInd belongs to edge
        ValidParentIndList = [];
        for j = 1:size(Node.ParentIndList,2)
            ParentInd = Node.ParentIndList(j);
            TrackEdgeList = cell2mat(NetObj.EdgeList(ParentInd, NodeInd));
            if (ismember(TrackInd, TrackEdgeList)~=0)
                ValidParentIndList = union(ValidParentIndList,ParentInd);
            end
        end 
        
        for i = 1:size(ValidParentIndList,2)
            ParentInd = ValidParentIndList(i);
            p_D_m1 = p_D(ParentInd,1);
            p_DT(TrackInd,NodeInd) = p_DT(TrackInd,NodeInd) + p_D_m1;
        end
    end
end
NetObj.p_DT = p_DT;
% Compute P_T matrix
p_T = ones(TrackNum, size(NetObj.NodeList,2));
p_T(:,1) = zeros(TrackNum,1);
for NodeInd=2:size(NetObj.NodeList,2)
    Node = NetObj.NodeList{NodeInd};
    for TrackInd = 1:TrackNum
        p_T(TrackInd, NodeInd) = p_U(NodeInd)*Li(Node.MeasInd, TrackInd)*p_DT(TrackInd,NodeInd);
    end
end
NetObj.p_T = p_T;
% Compute betta
betta = zeros(PointNum, TrackNum);
for MeasInd = 1:PointNum
    for TrackInd = 1:TrackNum
        % Get index list L_j of nodes in current measurement layer (MeasInd)
        L_j_Ind = find(cell2mat(cellfun(@(sas)sas.MeasInd, NetObj.NodeList, 'uni', false ))==MeasInd);
        for j = 1:size(L_j_Ind, 2)
            NodeInd = L_j_Ind(j);
            betta(MeasInd, TrackInd) = betta(MeasInd, TrackInd) + p_T(TrackInd,NodeInd);%p_U(NodeInd)*Li(MeasInd, TrackInd)*p_DT(TrackInd, NodeInd);
        end
    end
end

% Normalise
for j = 1:PointNum
    betta(j,:) = betta(j,:)/sum(betta(j,:),2);
end
betta = betta
betta_rescaled = betta./(betta(:,1)*ones(1,TrackNum));
betta_rescaled_reduced = betta_rescaled(:, 2:end);
betta_modified = [ones(1, TrackNum-1)*TrackNum-1; betta_rescaled_reduced];
betta_modified_transposed = betta_modified';
betta_transpose = betta_modified_transposed./(sum(betta_modified_transposed,2)*ones(1,PointNum+1))
%betta = betta
betta_trans = [ones(1,size(betta,2)-1)-sum(betta(:,2:end),1); betta(:,2:end)]'
NetObj.betta = betta;
NetObj.betta_trans = betta_transpose;
%end