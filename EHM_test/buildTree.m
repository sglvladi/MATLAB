%   buildTree.m                              Author: Lyudmil Vladimirov
%   ======================================================================>
%   Functionality: 
%       Build Track-Oriented Tree of association hypotheses and compute
%       the respective association probabilities (betta).
%       (normally executed for each cluster)
%
%   Input:
%       ValidationMatrix    - Matrix (m x n) containing all possible measurement
%                             to track associations (for cluster of interest).
%                             (Output from ObservationAssociation.m)
%       Li                  - Matrix (m x n) containing association likelihoods
%       (m: number of measurements, including dummy at index 1)
%       (n: number of tracks)

%   Output:
%       Structure NetObj:
%           #NetObj.NodeList - List of all Nodes (NodeObj) contained in net
%           #NetObj.EdgeList - Cell Matrix (n x n, n being the total number of
%                              nodes), where cell (i,j) contains the tracks
%                              contained in edge from parent node i, to
%                              child node j.
%           #NetObj.ValidationMatrix - Internal Validation matrix for
%                                      cluster.
%           #NetObj.p_D      - Computed p_D matrix
%           #NetObj.p_U      - Computed p_U matrix
%           #NetObj.p_DT     - Computed p_DT matrix
%           #NetObj.p_T      - Computed p_T matrix
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


function [TreeObj] = buildTree(ValidationMatrix, Li)
    
    % Get number of tracks/layers
    TrackNum = size(ValidationMatrix,1); 
    LayerNum = TrackNum+1; % Layer 1 is root layer

    % Augment ValidationMatrix
    %ValidationMatrix = [ones(1, TrackNum); ValidationMatrix']';

    % Get number of Maasurements
    PointNum = size(ValidationMatrix,2);

    % Define Node Object
    NodeObj.TrackInd = 0; % Track index
    NodeObj.MeasInd = []; % Measurement index
    NodeObj.ParentInd = []; % Parent index 
    NodeObj.ChildIndList = [];  % List of Child nodes
    NodeObj.Remainders = [1:PointNum]'; % Remaining measurements

    % Define Tree Object
    TreeObj.NodeList = [];
    sc = sum(ValidationMatrix,2);
    TreeObj.EdgeList = sparse(prod(sc), prod(sc));
    TreeObj.NodesVsLayer = 1; % n x 1 matrix, specifying the layer to which each node belongs 
                          %  (LayerNum = 1 is the dummy layer then virtually LayerNum = LayerNum+1)
    TreeObj.ValidationMatrix = ValidationMatrix;
    TreeObj.Li = Li;

    % Create Root Node
    TreeObj.NodeList{1} = NodeObj;

    % Build the tree
    % ==============>
    % For every track/layer
    for j=2:LayerNum
        
        % Get index list (L_jm1_Ind) of nodes in previous layer (j-1)
        L_jm1_Ind = find(TreeObj.NodesVsLayer == j-1); 

        % Get indeces of all associated measurements for track j
        M_j = [find(ValidationMatrix(j-1,:))]'; % i=0 for false alarm

        % For every node in L_jm1
        for i_jm1 = 1:size(L_jm1_Ind,1)

            % Index of parent node
            ParentInd = L_jm1_Ind(i_jm1);

            % Get all measurements to consider
            M_jm1 = union(1, intersect(M_j, TreeObj.NodeList{ParentInd}.Remainders)); 
            
            % For every measurement in M_jm1
            for i=1:size(M_jm1,1)

                % Get the measurement index
                MeasInd = M_jm1(i);
            
                % Create child node
                ChildInd = size(TreeObj.NodeList,2)+1;
                TreeObj.NodeList{ChildInd} = NodeObj;
                TreeObj.NodeList{ChildInd}.TrackInd = j-1;
                TreeObj.NodeList{ChildInd}.MeasInd = MeasInd;
                TreeObj.NodeList{ChildInd}.ParentInd = ParentInd;
                TreeObj.NodeList{ParentInd}.ChildIndList = union(TreeObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                TreeObj.NodesVsLayer(ChildInd,1) = j; % New node belongs to layer j

                % Create edge from parent to child
                TreeObj.EdgeList(ParentInd, ChildInd) = MeasInd;

                % Compute remainders
                TreeObj.NodeList{ChildInd}.Remainders = union(1,setdiff(TreeObj.NodeList{ParentInd}.Remainders, MeasInd)); 
            end
        end
    end
    
    % Compute association weights
    % ===========================>
    % Perform Backward-Forward algorithm
    %  *Reduced case of EHM BF, where we only have ONE parent per node 
    %   and ONE measurement per edge
    
    % Calculate the vector P_D
    p_D = zeros(size(TreeObj.NodeList,2),1);
    p_D(1,1) = 1;
    % for every node
    for NodeInd=2:size(TreeObj.NodeList,2)
        % Get Node object
        Node = TreeObj.NodeList{NodeInd};
        % Get Parent 
        ParentInd = Node.ParentInd;
        p_D_m1 = p_D(ParentInd,1);
        % Get measurement in Parent-Node edge
        MeasInd = TreeObj.EdgeList(ParentInd, NodeInd);
        % Calculate and store p_D for node
        p_D(NodeInd,1) =  Li(Node.TrackInd, MeasInd)*p_D_m1; 
    end
    TreeObj.p_D = p_D; 

    % Calculate the vector P_U
    p_U = zeros(size(TreeObj.NodeList,2),1);
    % Get index list (LeafIndList) of nodes in leaf layer and set p_U = 1
    LeafIndList = find(TreeObj.NodesVsLayer == LayerNum);
    p_U(LeafIndList,1) = 1;
    % for every Layer Starting from one to last (LayerNum-1)
    for LayerInd = LayerNum-1:-1:0
        % Get index list (L_j_Ind) of nodes in current layer
        L_j_Ind = find(TreeObj.NodesVsLayer == LayerInd);
        % for every node in the layer
        for i=1:size(L_j_Ind,1)
            p_U_temp = 0;
            % Get Node object
            NodeInd = L_j_Ind(i);
            Node = TreeObj.NodeList{NodeInd};
            % for every child of the node
            for j = 1:size(Node.ChildIndList,2)
                % Get Child
                ChildInd = Node.ChildIndList(j);
                ChildNode = TreeObj.NodeList{ChildInd};
                p_U_p1 = p_U(ChildInd,1);
                % Get measurement in Node-Child edge
                MeasInd = TreeObj.EdgeList(NodeInd, ChildInd);
                % Accumulate p_U_temp
                p_U_temp = p_U_temp + Li(ChildNode.TrackInd, MeasInd)*p_U_p1;
            end
            % Assign p_U
            p_U(NodeInd,1) = p_U_temp; 
        end
    end
    TreeObj.p_U = p_U;
    
    % Compute P_DT matrix
    p_DT = zeros(PointNum, size(TreeObj.NodeList,2));
    for NodeInd=1:size(TreeObj.NodeList,2)
        Node = TreeObj.NodeList{NodeInd};
        for MeasInd = 1:PointNum
            ParentInd = Node.ParentInd;
            if TreeObj.EdgeList(ParentInd, NodeInd)==MeasInd
                p_D_m1 = p_D(ParentInd,1);
                p_DT(MeasInd,NodeInd) = p_DT(MeasInd,NodeInd) + p_D_m1;
            end
        end
    end
    TreeObj.p_DT = p_DT;
    
    % Compute P_T matrix
    p_T = ones(PointNum, size(TreeObj.NodeList,2));
    p_T(:,1) = zeros(PointNum,1);
    for NodeInd=2:size(TreeObj.NodeList,2)
        Node = TreeObj.NodeList{NodeInd};
        for MeasInd = 1:PointNum
            p_T(MeasInd, NodeInd) = p_U(NodeInd)*Li(Node.TrackInd, MeasInd)*p_DT(MeasInd,NodeInd);
        end
    end
    TreeObj.p_T = p_T;

    % Compute betta
    betta = zeros(TrackNum, PointNum);
    for TrackInd = 1:TrackNum
        for MeasInd = 1:PointNum
            % Get index list L_j of nodes in current track layer (TrackInd+1)
            L_j_Ind = find( TreeObj.NodesVsLayer==TrackInd+1);
            for j = 1:size(L_j_Ind, 1)
                NodeInd = L_j_Ind(j);
                betta(TrackInd, MeasInd) = betta(TrackInd, MeasInd) + p_T(MeasInd,NodeInd);%p_U(NodeInd)*Li(TrackInd, TrackInd)*p_DT(TrackInd, NodeInd);
            end
        end
    end

    % Normalise
    for j = 1:TrackNum
        betta(j,:) = betta(j,:)/sum(betta(j,:),2);
    end
    TreeObj.betta = betta;
end