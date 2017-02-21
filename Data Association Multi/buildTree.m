function TreeObj = buildTree(parentNode, ClusterTracks, ValidationMatrix, TreeObj)
    firstTrackInd = ClusterTracks(1);
    GatedMeasurements = find(ValidationMatrix(:, firstTrackInd)');
    for i=1:size(GatedMeasurements,2)
        MeasInd = GatedMeasurements(i);
        %if(parentNode~=[0,0])
            %a = ismember([MeasInd,firstTrackInd], NetObj.NodeList{parentNode(1),parentNode(2)}.children);
            %if(a==0)
        %        NetObj.NodeList{parentNode(1),parentNode(2)}.children = [NetObj.NodeList{parentNode(1),parentNode(2)}.children; [MeasInd,firstTrackInd]];;
        
            %end
        %end
        BranchInd = size(TreeObj.NodeList{MeasInd,firstTrackInd}.branches,2)+1; 
        TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd} = [parentNode];
        if(size(ClusterTracks,2)>1)
            nextTrackInd = ClusterTracks(2);
            %NetObj.NodeList{MeasInd,firstTrackInd}.children = [NetObj.NodeList{MeasInd,firstTrackInd}.children; [MeasInd,nextTrackInd]];
            nextValidationMatrix = ValidationMatrix;
            nextValidationMatrix(MeasInd,:) = 0;
            nextGatedMeasurements = find(nextValidationMatrix(:, nextTrackInd)');
            for i=1:size(nextGatedMeasurements,2)
                nextMeasInd = nextGatedMeasurements(i);
                TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd} = [TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd}; [nextMeasInd, nextTrackInd]];
            end
            TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd} = [TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd}; [5, nextTrackInd]];
            TreeObj = buildTree([MeasInd, firstTrackInd], ClusterTracks(ClusterTracks~=firstTrackInd), nextValidationMatrix, TreeObj);
        else
            TreeObj.LeafNodes = [TreeObj.LeafNodes; [MeasInd, firstTrackInd]];
        end
        
    end
    MeasInd = size(ValidationMatrix,1)+1;
    %if(parentNode~=[0,0])
        %a = ismember([MeasInd,firstTrackInd], NetObj.NodeList{parentNode(1),parentNode(2)}.children);
        %if(a==0)
          %  NetObj.NodeList{parentNode(1),parentNode(2)}.children = [NetObj.NodeList{parentNode(1),parentNode(2)}.children; [MeasInd,firstTrackInd]];;
    
        %end
    %end
    BranchInd = size(TreeObj.NodeList{MeasInd,firstTrackInd}.branches,2) +1; 
    TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd} = [parentNode];
    if(size(ClusterTracks,2)>1)
        nextTrackInd = ClusterTracks(2);
        %NetObj.NodeList{MeasInd,firstTrackInd}.children = [NetObj.NodeList{MeasInd,firstTrackInd}.children; [MeasInd,nextTrackInd]];
        nextValidationMatrix = ValidationMatrix;
        nextGatedMeasurements = find(nextValidationMatrix(:, nextTrackInd)');
        for i=1:size(nextGatedMeasurements,2)
            nextMeasInd = nextGatedMeasurements(i);
            TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd} = [TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd}; [nextMeasInd, nextTrackInd]];
        end
        TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd} = [TreeObj.NodeList{MeasInd,firstTrackInd}.branches{BranchInd}; [5, nextTrackInd]];
        TreeObj = buildTree([MeasInd, firstTrackInd], ClusterTracks(ClusterTracks~=firstTrackInd), nextValidationMatrix, TreeObj);
    else
        TreeObj.LeafNodes = [TreeObj.LeafNodes; [MeasInd, firstTrackInd]];
    end
    
end