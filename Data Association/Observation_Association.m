function [TrackList, Validation_matrix] = Observation_Association(TrackList,DataList,EKF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Association - connects between data and tracks 
% Input:
%   TrackList    - kalman structure list and more
% Output:
%   TrackList    - kalman structure list updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
GateLevel   = 5; % 98.9% of data in gate
TrackNum    = size(TrackList,2);
PointNum    = size(DataList,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check un-initalized tracks (state 0)
%for i=1:TrackNum,
%    if TrackList{i}.TrackObj.State == Par.State_Undefined, error('Undefined Tracking object'); end;
%end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find associated data
% init association matrix
DistM = ones(TrackNum,PointNum)*1000;
for i=1:TrackNum,
    
    % extract track coordinates and covariance
    %x = TrackList{i}.TrackObj.x;
    %P = TrackList{i}.TrackObj.P;
    %F = TrackList{i}.TrackObj.F;
    %H = TrackList{i}.TrackObj.H;
    %Q = TrackList{i}.TrackObj.Q;
    %R = TrackList{i}.TrackObj.R;
    TrackList{i}.TrackObj = EKF.Predict(TrackList{i}.TrackObj);
    % make prediction
    %xpred   = F*x;
    %Ppred   = F*P*F' + Q;
    %ypred   = H*xpred;
    %S       = H*Ppred*H' + R;
    %x_pred
    % measure points
    for j=1:PointNum,

        % distance
        DistM(i,j)  = mahalDist(DataList(:,j), TrackList{i}.TrackObj.z_pred, TrackList{i}.TrackObj.S, 2); 
        
    end;
    
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thresholding/ gating
Validation_matrix = DistM < GateLevel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% association
%for i=1:TrackNum,
    
    % if is empty then DataInd is empty
%    ValidAssociatedInd              = find(DistLabels(i,:))
%    TrackList{i}.TrackObj.DataInd   = ValidAssociatedInd;
 
%end;




