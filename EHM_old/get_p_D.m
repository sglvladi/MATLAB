function p_D = get_p_D(NodeInd, NetObj, Li)
    Node = NetObj.NodeList{NodeInd}; % current Node
    p_D = 0; % Initiate p_D
    
    if(~isempty(Node.ParentIndList))
        % for all Parents
        for j = 1:size(Node.ParentIndList,2)
            p_Dtemp = 0; % Initiate temporary p_D
            ParentInd = Node.ParentIndList(j); % Index of currently considered Parent
            ParentNode = NetObj.NodeList{NodeInd}; % current Parent node
            MeasInd = ParentNode.MeasInd; % Measurement/Layer Index
            TrackEdgeList = cell2mat(NetObj.EdgeList(ParentInd,NodeInd)); %List of possibly associated track indices

            % for each track in the list of possibly associated tracks
            % ========================================================
            % => Possibly replace with:
            %
            %       p_Dtemp = size(TrackEdgeList,2) * Li(MeasInd,TrackInd) * get_p_D(ParentInd, NetObj, Li);
            %
            %    , since the path followed to root is the same for each track
            %    inside the edge.
            for i = 1:size(TrackEdgeList,2)
                TrackInd = TrackEdgeList(i); % Get the track index
                p_Dtemp = p_Dtemp + Li(MeasInd,TrackInd) * get_p_D(ParentInd, NetObj, Li);
            end

            % Sum for all parent nodes
            p_D = p_D + p_Dtemp;
        end
    else
        p_D = 1;
    end
end