function TrackList = Structure_PDAF_Init_Track(Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Init_Track - initializes the Kalman filter
% matrix with respect to dimensions and other track related data
% Input:
%  Par           - parameters for to intialization
% Output:
%   TrackList    - kalman structure list and more
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default
% all the parameters must be defined in Par structure

% tracker structure
[F,H,Q,R,ObservInd,initx,initP] = Kalman_Filter_Init(Par.dT,Par.ModelDim,Par.ProblemDim,Par.StateVar,Par.ObserVar);
TrackObj.F          = F;
TrackObj.H          = H;
TrackObj.Q          = Q;
TrackObj.R          = R;
TrackObj.ObservInd  = ObservInd;

% initial location and variance
TrackObj.x          = initx;
TrackObj.P          = initP;

% intialization state
TrackObj.State      = Par.State_FisrtInit;  
% life time counter
TrackObj.LifeTime   = 0;  
% log likelihood
TrackObj.LogLike    = Par.GateLevel;  % is not zero but it is not so bad start
% history
TrackObj.Hist       = zeros(size(F,1),Par.HistLen);

% distribute the trackers evently over the measurement space
dy1                 = diff(Par.Y1Bounds);
dy2                 = diff(Par.Y2Bounds);
dy3                 = diff(Par.Y3Bounds);
switch Par.ProblemDim,
     case 2,

        TrackNumY1          = max([1 floor(sqrt(dy1/dy2*Par.TrackNum))]);
        TrackNumY2          = ceil(Par.TrackNum/TrackNumY1);
        % generates grid of centers
        [yy1,yy2]           = meshgrid(((Par.Y1Bounds(1)+dy1/TrackNumY1/2):dy1/TrackNumY1:Par.Y1Bounds(2)),...
                                       ((Par.Y2Bounds(1)+dy2/TrackNumY2/2):dy2/TrackNumY2:Par.Y2Bounds(2)));
        CenterData          = [yy1(:) yy2(:)]';
    case 3,
        gridNum             = ceil(Par.TrackNum^(1/3));
        % generates grid of centers
        [yy1,yy2,yy3]       = meshgrid( linspace(Par.Y1Bounds(1),Par.Y1Bounds(2),gridNum),...
                                        linspace(Par.Y2Bounds(1),Par.Y1Bounds(2),gridNum),...
                                        linspace(Par.Y3Bounds(1),Par.Y1Bounds(2),gridNum));
        CenterData          = [yy1(:) yy2(:) yy3(:)]';
        
        
end

% init the list
TrackList = [];
for i=1:Par.TrackNum,
    
    TrackList{i}.TrackObj = TrackObj;
    TrackList{i}.TrackObj.x(ObservInd) = CenterData(:,i);

end;


