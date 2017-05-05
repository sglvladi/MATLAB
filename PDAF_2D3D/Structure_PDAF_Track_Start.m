function TrackList = Structure_PDAF_Track_Start(TrackList,DataList,Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Track_Start - initialization of new tracks
% Input:
%   TrackList    - kalman structure list and more
%   DataList     - 2 x DataPointsNum contains the relevant data for time 1
%   Par          - parameters
% Output:
%   TrackList    - kalman structure list and more
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
TrackNum    = Par.TrackNum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find only the valid data  (no NaN)
UnAssociatedData = sum(isnan(DataList)) == 0;       % 1 - is valid data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find an unassociated data (No track is using)
% label an undefined tracks
UnDefinedTrackLabel = zeros(1,TrackNum);

% loop throught all tracks
for i=1:TrackNum,
    
    % only valid tracks - in Track and Init states
    ValidStates = [Par.State_Track Par.State_FisrtInit:Par.State_LastInit];
    if any(TrackList{i}.TrackObj.State == ValidStates), 
        UnAssociatedData(TrackList{i}.TrackObj.DataInd) = 0; % 1 - is Unassociated
    end;
    
    % only undefined states
    if TrackList{i}.TrackObj.State == Par.State_Undefined, 
        UnDefinedTrackLabel(i) = 1; % 1 - is undefined
    end;
    
end;

UnAssocDataInd = find(UnAssociatedData);
UnAssocDataLen = length(UnAssocDataInd);
if Par.ShowOn, fprintf('Unassociated points number: %d\n',UnAssocDataLen); end;

UnDefinedTrackInd = find(UnDefinedTrackLabel);
UnDefinedTrackLen = length(UnDefinedTrackInd);
if Par.ShowOn, fprintf('Undefined state tracks number: %d\n',UnDefinedTrackLen); end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare the starting locations for each track 
% and resolve different situations
if UnDefinedTrackLen == 0,      % no free tracks
    
    if Par.ShowOn,  disp('No undefined tracks'); end;
    
    return;
end;

% there are more tracks then data
if UnAssocDataLen < UnDefinedTrackLen, % there is more tracks that unassociated data points
    % !!! SIMPLE RANDOM INIT - can be improved
switch Par.ProblemDim,
     case 2,
        dy1 = diff(Par.Y1Bounds);
        dy2 = diff(Par.Y2Bounds);
        RandLocY1 = rand(1,UnDefinedTrackLen-UnAssocDataLen).*dy1 + Par.Y1Bounds(1);
        RandLocY2 = rand(1,UnDefinedTrackLen-UnAssocDataLen).*dy2 + Par.Y2Bounds(1);
        RandLoc   = [RandLocY1;RandLocY2];
        DataLoc   = [RandLoc DataList(:,UnAssocDataInd)];
     case 3,
        dy1 = diff(Par.Y1Bounds);
        dy2 = diff(Par.Y2Bounds);
        dy3 = diff(Par.Y3Bounds);
        RandLocY1 = rand(1,UnDefinedTrackLen-UnAssocDataLen).*dy1 + Par.Y1Bounds(1);
        RandLocY2 = rand(1,UnDefinedTrackLen-UnAssocDataLen).*dy2 + Par.Y2Bounds(1);
        RandLocY3 = rand(1,UnDefinedTrackLen-UnAssocDataLen).*dy3 + Par.Y3Bounds(1);
        RandLoc   = [RandLocY1;RandLocY2;RandLocY3];
        DataLoc   = [RandLoc DataList(:,UnAssocDataInd)];
 end    
else
    % !!! check if the index permutation is usefull
    DataLoc   = DataList(:,UnAssocDataInd(1:UnDefinedTrackLen));    
end;

%fprintf(' %d tracks are initiated and %d from random location\n',UnDefinedTrackLen,UnAssocDataLen);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop throught all undefined tracks
for i=1:UnDefinedTrackLen,
         
    % reset
    xnew    = zeros(size(TrackList{UnDefinedTrackInd(i)}.TrackObj.x));
    xnew(TrackList{UnDefinedTrackInd(i)}.TrackObj.ObservInd) = DataLoc(:,i);
    Pnew    = TrackList{UnDefinedTrackInd(i)}.TrackObj.Q;
    
    % Reset different Kalman structures
    TrackList{UnDefinedTrackInd(i)}.TrackObj.x     = xnew;
    TrackList{UnDefinedTrackInd(i)}.TrackObj.P     = Pnew * 10;
    TrackList{UnDefinedTrackInd(i)}.TrackObj.State = Par.State_FisrtInit;
    TrackList{UnDefinedTrackInd(i)}.TrackObj.LifeTime = 0;
    TrackList{UnDefinedTrackInd(i)}.TrackObj.LogLike  = Par.GateLevel;
    
    % history reset
    TrackList{UnDefinedTrackInd(i)}.TrackObj.Hist  = zeros(size(TrackList{UnDefinedTrackInd(i)}.TrackObj.Hist));
    
end;
