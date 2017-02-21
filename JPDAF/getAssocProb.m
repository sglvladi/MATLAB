function betta = getAssocProb(MeasInd, TrackInd, HypMap, HypProb)
    ValidHyps = [];
    keys(HypMap);
    size(keys(HypMap));
    
    % get the hypotheses where the node [MeasInd, TrackInd] exists
    for i = 1:size(keys(HypMap),2)-1
        a = ismember([MeasInd, TrackInd], HypMap(i), 'rows');
        if (a==1)
            ValidHyps = [ValidHyps, i]; % Add to Valid hypotheses
        end
    end
    
    betta = 0;
    % sum the Probabilities for all Valid hypothses
    for i = 1:size(ValidHyps,2)
        betta = betta + HypProb(ValidHyps(i));
    end
    
    % Normalise betta
    betta = betta/(sum(HypProb(1:end-1)));
    
end
