function [HypMap, HypProb] = buildHyp(HypMap, HypProb, Li, ClusterTracks, ValidationMatrix, branch, bettaNTFA)
    firstTrackInd = ClusterTracks(1);
    %HypMap = containers.Map;
    GatedMeasurements = find(ValidationMatrix(:, firstTrackInd)');
    
    %For all associated measurements (except FA)
    for i=1:size(GatedMeasurements,2)
        branchTemp = branch; % temporary branch
        MeasInd = GatedMeasurements(i); % Index of current measurement
        Hyps = keys(HypMap); % list of hypothesis access keys
        lastHyp = size(Hyps,2); % Number of different hypothesises
        
        % Extend temporary branch and update Hypothesis probability
        %  for new (extended) branch.
        branchTemp = [branchTemp; [MeasInd, firstTrackInd]];
        HypProb(lastHyp) = HypProb(lastHyp)*Li(MeasInd, firstTrackInd);
        
        % If there exist other unconsidered tracks/members of the cluster
        %  then keep extending branch
        if(size(ClusterTracks,1)>1)
            nextTrackInd = ClusterTracks(2);
            %NetObj.NodeList{MeasInd,firstTrackInd}.children = [NetObj.NodeList{MeasInd,firstTrackInd}.children; [MeasInd,nextTrackInd]];
            nextValidationMatrix = ValidationMatrix;
            nextValidationMatrix(MeasInd,:) = 0;
            nextGatedMeasurements = find(nextValidationMatrix(:, nextTrackInd)');
            [HypMap, HypProb] = buildHyp(HypMap, HypProb, Li, ClusterTracks(ClusterTracks~=firstTrackInd), nextValidationMatrix, branchTemp, bettaNTFA);
        % else, stop extending, store branch as a complete hypothesis
        %  end create new entry for next hypothesis
        % NOTE: Creating the new entry only when the previous hypothesis
        %       has been  
        else
            HypMap(lastHyp) = branchTemp;
            HypMap(lastHyp+1) = [];
            HypProb(lastHyp+1) = bettaNTFA;
            Hyps = keys(HypMap);
            lastHyp = size(Hyps,2);
        end
    end
    
    % Add FA branch
    MeasInd = size(ValidationMatrix,1)+1;
    Hyps = keys(HypMap);
    lastHyp = size(Hyps,2);
    branch = [branch; [MeasInd, firstTrackInd]];
    HypProb(lastHyp) = HypProb(lastHyp)*Li(MeasInd, firstTrackInd);
    if(size(ClusterTracks,1)>1)
        nextTrackInd = ClusterTracks(2);
        %NetObj.NodeList{MeasInd,firstTrackInd}.children = [NetObj.NodeList{MeasInd,firstTrackInd}.children; [MeasInd,nextTrackInd]];
        nextValidationMatrix = ValidationMatrix;
        nextGatedMeasurements = find(nextValidationMatrix(:, nextTrackInd)');
        [HypMap, HypProb] = buildHyp(HypMap, HypProb, Li, ClusterTracks(ClusterTracks~=firstTrackInd), nextValidationMatrix, branch, bettaNTFA);
    else
        HypMap(lastHyp) = branch;
        HypMap(lastHyp+1) = [];
        HypProb(lastHyp+1) = bettaNTFA;
        Hyps = keys(HypMap);
        lastHyp = size(Hyps,2);
    end
    
end