%   Observation_Association.m                   Author: Lyudmil Vladimirov
%   ======================================================================>
%   Functionality: 
%       Perform gate validation to compute validation matrix and NT/FA
%       density.
%   
%   Input: 
%       TrackList    - List of all target tracks at time k(TrackObj's)
%       DataList     - List of all measurements at time k
%       
%   Output:
%       ValidationMatrix    - Matrix containing all possible measurement to
%                             track associations.
%       bettaNTFA           - New track/False alarm density (assumed to be same)
%       
%   
%   Dependencies: mahalDist.m 
%   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [TrackList, ValidationMatrix, bettaNTFA] = Observation_Association(TrackList,DataList,UKF)


% parameters
GateLevel   = 5; % 98.9% of data in gate
TrackNum    = size(TrackList,2);
PointNum    = size(DataList,2);
ObsDim      = size(DataList,1);
C = pi;     % volume of the 2-dimensional unit hypersphere (change for different Dim no)


% variables
tot_gate_area = 0;

% Initiate Mahalanobis distance storage matrix
DistM = ones(TrackNum,PointNum)*1000;

% for each track
for i=1:TrackNum
    
    % Perform state and measurement prediction
    TrackList{i}.TrackObj = UKF.Predict(TrackList{i}.TrackObj);
    
    % Compute total gate area
    tot_gate_area = tot_gate_area + C*GateLevel^(ObsDim/2)*det( TrackList{i}.TrackObj.S)^(1/2);
    
    %for each measurement
    for j=1:PointNum

        % Compute and store mahalanobis distance
        DistM(i,j)  = mahalDist(DataList(:,j), TrackList{i}.TrackObj.z_pred, TrackList{i}.TrackObj.S, 2); 
        
    end;
    
end;

% Apply gate threshold to obtain the Validation Matrix
ValidationMatrix = DistM < GateLevel;

% Compute New Track/False Alarm density
bettaNTFA = sum(ValidationMatrix(:))/tot_gate_area;

end





