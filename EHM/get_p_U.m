function p_U = get_p_U(NodeInd, NetObj, Li)
    Node = NetObj.NodeList{NodeInd}; % current Node
    p_U = 0; % Initiate p_D
    
    if(~isempty(Node.ChildIndList))
        % for all Children
        for j = 1:size(Node.ChildIndList,2)
            p_Utemp = 0; % Initiate temporary p_U
            ChildInd = Node.ChildIndList(j); % Index of currently considered Parent
            ChildNode = NetObj.NodeList{ChildInd}; % current Parent node
            MeasInd = ChildNode.MeasInd; % Measurement/Layer Index
            TrackEdgeList = cell2mat(NetObj.EdgeList(NodeInd,ChildInd)); %List of possibly associated track indices

            % for each track in the list of possibly associated tracks
            for i = 1:size(TrackEdgeList,2)
                TrackInd = TrackEdgeList(i); % Get the track index
                p_Utemp = p_Utemp + Li(MeasInd,TrackInd) * get_p_U(ChildInd, NetObj, Li);
            end
            p_U = p_U + p_Utemp;
        end
    else
        p_U = 1;
    end
end