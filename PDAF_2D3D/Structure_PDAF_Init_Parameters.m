function Par = Structure_PDAF_Init_Parameters();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Init_Parameters - initializes different parameters of 
% the algorithm
% Input:
%   -
% Output:
%   Par          - parameter structure 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% general management
Par.ShowOn          = 0;        % show different information
Par.RecordOn        = 0;        % record the movie


% init data and tracks
Par.TrajIndex       = [12 11 8 8 8 1]; %ones(1,4).*8; %[1 5 7 8]; %ones(1:2).*8; %[1 2];   % max 7
Par.PointNum        = 18;       % number of clutter and data points must be greater then TrajIndex length
Par.NaNDensity      = 0.05;        %.1; 
Par.Y1Bounds        = [0 1];    % approximate bounds for measurements in Y1
Par.Y2Bounds        = [0 1];    % approximate bounds for measurements in Y2
Par.Y3Bounds        = [0 1];    % approximate bounds for measurements in Y2
Par.TrajNum         = length(Par.TrajIndex); % save this for future use
Par.Nv              = 0.001;     % noise variance
Par.dT              = 1/30;     % time between measurements in seconds
Par.Time            = 5;        % simulation stopping time in seconds

% Kalman filter properties
Par.StateVar        = (0.5)^2;  % state variance sqrt(StateVar)/dT
Par.ObserVar        = (0.01)^2;  % observation variance sqrt(StateVar)/dT

% Track properties
Par.TrackNum        = 12;        % number of trackers
Par.ProblemDim      = 2;        % ProbDim     - dimensionality of the problem (x, x-y or x-y-z)
Par.ModelDim        = 2;        % ModelDim    - constant velocity or constant acceleration
Par.HistLen         = 3;        % HistLen     - number of past state history for each tracker
Par.HistGateLevel   = 0.1;      % determines history separation value
Par.LogLikeAlpha    = 0.3;      % log likelihood forget factor

%Par.GateLev        = 9.2;      % - 99 percent in 2D
Par.GateLevel       = 5;        % - 98.9 percent
%Par.GateLevel       = 4;       % - 86.5 percent

% PDAF parameters
Par.UsePDAF         = 1 ;       % switch to the PDAF mode
Par.PG              = .9;       % probabiliry of gating
Par.PD              = .8;       % probability of detection;

% Tracker states
Par.State_Undefined = 0;        % track is a free resources and not defined
Par.State_FisrtInit = 1;        % track is in the first initalization state
Par.State_LastInit  = 3;        % track is in the last initalization state
Par.State_Track     = 10;       % track is in the tracking mode
Par.State_FirstCoast= 20;       % track is in first cost state
Par.State_LastCoast = 25;       % track is in last cost state

