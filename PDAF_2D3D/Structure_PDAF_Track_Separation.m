function TrackList = Structure_PDAF_Track_Separation(TrackList,DataList,Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Track_Separation - identifies tracks with the
% same track history and discards them.
% Input:
%   TrackList    - kalman structure list and more
%   DataList     - 2 x DataPointsNum contains the relevant data for time (is not used here)
%   Par          - parameters
% Output:
%   TrackList    - kalman structure list and more
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
TrackNum    = Par.TrackNum;
% proportional to the size of the tracking region
HistGate    = Par.HistGateLevel*diff(Par.Y1Bounds)*diff(Par.Y2Bounds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine the distance using history data
HistDist = ones(TrackNum,TrackNum)*10000;
for i=1:TrackNum-1, 
    for j=i+1:TrackNum, 
        HistDist(i,j) = sum(std(TrackList{i}.TrackObj.Hist - TrackList{j}.TrackObj.Hist)); 
    end; 
end;

SameTracks = HistDist < HistGate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple decision about track separation:
% 1.Find the closest tracks and dump them.
% It is improved using additional logic because it preferes
% the tracks from top of the list

% determine only the valid tracks and their LifeTimes
ValidTrackLabel = zeros(1,TrackNum);
TracksLifeTime  = zeros(1,TrackNum);

% loop throught all tracks
for i=1:TrackNum,
    
    % only valid tracks - not in Undefined state
    if TrackList{i}.TrackObj.State ~= Par.State_Undefined, 
        ValidTrackLabel(i) = 1; % 1 - is OK
        
        % store life time
        TracksLifeTime(i) = TrackList{i}.TrackObj.LifeTime;
    end;  
end;

ValidTrackInd = find(ValidTrackLabel);
ValidTrackLen = length(ValidTrackInd);
if ValidTrackLen == 0, 
    if Par.ShowOn, disp('Undefined States - check!!!'); end;
    return;
end;
if Par.ShowOn, fprintf('Valid state tracks number: %d\n',ValidTrackLen); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort valid tracks by their life times
[unused_values,sorted_ind] = sort(TracksLifeTime(ValidTrackInd));
SortedTrackInd = ValidTrackInd(sorted_ind(end:-1:1)); % descending order

% chnge matrix rows and discard Undefined tracks
SameTracks = SameTracks(SortedTrackInd,:);


% loop throught valid tracks
for i=1:ValidTrackLen,
    
    % determine the closest tracks to track i
    SameTrackInd = find(SameTracks(i,:));

    for j=1:length(SameTrackInd), % if isempty(SameTrackInd), wont be executed
        
        % goes to undefined state
        TrackList{SameTrackInd(j)}.TrackObj.State = Par.State_Undefined;
        % reset the row of distance matrix
        SameTracks(SameTrackInd(j),:) = 0;
    end;
end;    
    

 