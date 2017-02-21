% Create ValidationMatrix
ValidationMatrix = [1 1 0 1; 1 0 1 0; 1 1 0 0;  1 0 0 1; 1 1 1 0];

% Get number of Track and measurements
TrackNum = size(ValidationMatrix,2);
PointNum = size(ValidationMatrix,1); %j=1 is a dummy measurement/layer
LayerNum = TrackNum; % Layer 1 is root layer

% Define Node Object
NodeObj.MeasInd = 0; % Measurement index
NodeObj.TrackIndList = []; % List of associated tracks
NodeObj.ParentIndList = []; % List of Parent nodes
NodeObj.ChildIndList = [];  % List of Child nodes
NodeObj.Remainders = [ 1:TrackNum ]'; % Remaining tracks

% Define Net Object
NetObj.NodeList = [];
NetObj.EdgeList = [];

% Define Root Node
NetObj.NodeList{1} = NodeObj;


% For every measurement/layer
for j=1:LayerNum
    
    % Get index list (L_jm1_Ind) of nodes in previous layer (j-1)
    L_jm1_Ind = find(cell2mat(cellfun(@(sas)sas.MeasInd, NetObj.NodeList, 'uni', false))==j-1);
    %L_jm1 = NetObj.NodeList{j-1,:};
    
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
            
            % Get index list L_j of nodes in current layer (j)
            %L_j = NetObj.NodeList{j,:};
            L_j_Ind = find(cell2mat(cellfun(@(sas)sas.MeasInd, NetObj.NodeList, 'uni', false ))==j);
            
            % if current Layer is empty
            if isempty(L_j_Ind)
                % Index of child node
                ChildInd = size(NetObj.NodeList,2)+1;
                
                % Create child node
                NetObj.NodeList{ChildInd} = NodeObj;
                NetObj.NodeList{ChildInd}.MeasInd = j;
                NetObj.NodeList{ChildInd}.TrackIndList = union(NetObj.NodeList{L_jm1_Ind(i_jm1)}.TrackIndList(:,:), T_jm1(i));
                NetObj.NodeList{ChildInd}.TrackIndList = NetObj.NodeList{ChildInd}.TrackIndList(:);
                NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);
               
                % Create edge from parent to child
                NetObj.EdgeList{ParentInd, ChildInd} =[TrackInd];
                
                % Compute remainders
                T_rem_j = []; % T_j+1:mk[N^(j)_(i_j)]
                for j_sub = j+1:LayerNum
                    T_rem_j = union(T_rem_j, find(ValidationMatrix(j_sub,:))');
                end
                NetObj.NodeList{ChildInd}.Remainders = setdiff(T_rem_j,setdiff(getAncestors(ChildInd,[],NetObj.EdgeList),1)); 
            else
                %Possibly replace with R_jm1 = setdiff(NetObj.NodeList{ParentInd}.Remainders, TrackInd);
                T_rem_jm1_ti = []; % T_j+1:mk[N^(j-1)_(i_j-1),t_i]
                for j_sub = j:LayerNum
                    T_rem_jm1_ti = union(T_rem_jm1_ti, find(ValidationMatrix(j_sub,:))');
                end
                
                % Possibly replace with R_jm1 = setdiff(NetObj.NodeList{ParentInd}.Remainders, TrackInd);
                %R_jm1 = setdiff(T_rem_jm1_ti,setdiff(union(union(getAncestors(ParentInd,[],NetObj.EdgeList),NetObj.NodeList{ParentInd}.TrackIndList), TrackInd),1));
                R_jm1 = union(setdiff(NetObj.NodeList{ParentInd}.Remainders, TrackInd),1);
                R_jm1 = R_jm1(:); % Enforce that R_jm1 is a column vector
               
                % For all nodes in current layer
                for i_j=1:size(L_j_Ind,2)
                    ChildInd = L_j_Ind(i_j);
                    
                    % Compute remainders
                    T_rem_j = []; % T_j+1:mk[N^(j)_(i_j)]
                    for j_sub = j+1:LayerNum
                        T_rem_j = union(T_rem_j, find(ValidationMatrix(j_sub,:))');
                    end
                    %R_j = setdiff(T_rem_j,setdiff(union(getAncestors(ChildInd,[],NetObj.EdgeList),NetObj.NodeList{ChildInd}.TrackIndList),1)); 
                    R_j = NetObj.NodeList{ChildInd}.Remainders;
                     
                    % If the list of node's list of remainders is equal to R_jm1
                    if (isequal(R_jm1,R_j)||j==LayerNum) 
                        NetObj.NodeList{ChildInd}.TrackIndList = union(NetObj.NodeList{ChildInd}.TrackIndList,TrackInd);
                        if size(NetObj.EdgeList,1)<ParentInd
                            NetObj.EdgeList{ParentInd, ChildInd} = TrackInd;
                        else
                            NetObj.EdgeList{ParentInd, ChildInd} = union(NetObj.EdgeList{ParentInd, ChildInd}, TrackInd);
                        end
                        NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                        NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);
               
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
                    NetObj.NodeList{ChildInd}.TrackIndList = NetObj.NodeList{ChildInd}.TrackIndList(:);
                    NetObj.NodeList{ParentInd}.ChildIndList = union(NetObj.NodeList{ParentInd}.ChildIndList, ChildInd);
                    NetObj.NodeList{ChildInd}.ParentIndList = union(NetObj.NodeList{ChildInd}.ParentIndList, ParentInd);
               
                    % Create edge from parent to child
                    NetObj.EdgeList{ParentInd, ChildInd} =[TrackInd];
                    
                    % Compute remainders
                    T_rem_j = []; % T_j+1:mk[N^(j)_(i_j)]
                    for j_sub = j+1:LayerNum
                        T_rem_j = union(T_rem_j, find(ValidationMatrix(j_sub,:))');
                    end
                    R_j = setdiff(T_rem_j,setdiff(getAncestors(ChildInd,[],NetObj.EdgeList),1));
                    NetObj.NodeList{ChildInd}.Remainders = R_j; 
                end
            end

        end
    end
end    