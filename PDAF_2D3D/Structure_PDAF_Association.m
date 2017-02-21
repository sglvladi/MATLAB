function TrackList = Structure_PDAF_Association(TrackList,DataList,Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Association - connects between data and tracks 
% Input:
%   TrackList    - kalman structure list and more
%   DataList     - 2 x DataPointsNum contains the relevant data for time t
%   Par          - parameters
% Output:
%   TrackList    - kalman structure list updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
GateLevel   = Par.GateLevel^Par.ProblemDim;
TrackNum    = Par.TrackNum;
PointNum    = Par.PointNum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check un-initalized tracks (state 0)
for i=1:TrackNum,
    if TrackList{i}.TrackObj.State == Par.State_Undefined, error('Undefined Tracking object'); end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find missing data
ValidDataLabel = zeros(1,PointNum);
for i=1:PointNum,
    
    % not nan
    ValidDataLabel(i)= all(~isnan(DataList(:,i)))

end;
ValidDataInd = find(ValidDataLabel);
ValidDataNum = length(ValidDataInd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resolve zeros
if ValidDataNum == 0, 
    if Par.ShowOn, 
        fprintf('No valid data.\n'); % 
    end;
    ResolvedValidDataNum = 1;   % to prevent zero columns DistM matrix
else
    ResolvedValidDataNum = ValidDataNum;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find associated data
% init association matrix
DistM = ones(TrackNum,ResolvedValidDataNum)*1000;
for i=1:TrackNum,
    
    % extract track coordinates and covariance
    x = TrackList{i}.TrackObj.x;
    P = TrackList{i}.TrackObj.P;
    F = TrackList{i}.TrackObj.F;
    H = TrackList{i}.TrackObj.H;
    Q = TrackList{i}.TrackObj.Q;
    R = TrackList{i}.TrackObj.R;

    % make prediction
    xpred   = F*x;
    Ppred   = F*P*F' + Q;
    ypred   = H*xpred;
    S       = H*Ppred*H' + R;
    
    % measure points
    % if ValidDataNum=0 - skipped
    for j=1:ValidDataNum,
        
        % extract measured points
        ym = DataList(:,ValidDataInd(j));

        % distance
        DistM(i,j)  = gaussian_prob(ym, ypred, S, 2); 
        
    end;
    
end;
DistM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thresholding/ gating
DistLabels = DistM < GateLevel
%DistM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% association
for i=1:TrackNum,
    
    % if is empty then DataInd is empty
    ValidAssociatedInd              = find(DistLabels(i,:));
    TrackList{i}.TrackObj.DataInd   = ValidDataInd(ValidAssociatedInd);
 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unassociated data and new tracks
UnAssocDataInd = find(sum(DistLabels)==0);
UnAssocDataLen = length(UnAssocDataInd);
if Par.ShowOn, fprintf('Unassociated points number: %d\n',UnAssocDataLen); end;
UnAssocTrackInd = find(sum(DistLabels,2)==0);
UnAssocTrackLen = length(UnAssocTrackInd);
if Par.ShowOn, fprintf('Unassociated tracks number: %d\n',UnAssocTrackLen); end;




