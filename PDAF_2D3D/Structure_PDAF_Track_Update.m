function TrackList = Structure_PDAF_Track_Update(TrackList,DataList,Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Track_Update - performs track update
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
alpha       = Par.LogLikeAlpha;
PG          = Par.PG;
PD          = Par.PD;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop throught all tracks
for i=1:TrackNum,
      
    %------------------------------------------------
    % extract track coordinates, covariance and data
    x       = TrackList{i}.TrackObj.x;
    P       = TrackList{i}.TrackObj.P;
    F       = TrackList{i}.TrackObj.F;
    H       = TrackList{i}.TrackObj.H;
    Q       = TrackList{i}.TrackObj.Q;
    R       = TrackList{i}.TrackObj.R;
    DataInd = TrackList{i}.TrackObj.DataInd;
    State   = TrackList{i}.TrackObj.State;
    Hist    = TrackList{i}.TrackObj.Hist;
    
    %------------------------------------------------
    % associated data check
    if isempty(DataInd), 
        
        % !!! SWITCH didn't work
        if State == Par.State_Undefined,
            error('Undefined track');
        elseif any(State == (Par.State_FisrtInit:Par.State_LastInit)),      % initialization states
            NextState = Par.State_Undefined;                % reset the initialization 
        elseif any(State == (Par.State_FirstCoast:Par.State_LastCoast-1)),  % coast mode states
            NextState = State + 1;                          % naxt coast mode states
        elseif State == Par.State_LastCoast,                % final coast mode state
            NextState = Par.State_Undefined;                % reset the track 
        else                                                % track state
            NextState = Par.State_FirstCoast;               % first coast mode state
        end;
        
    else      
        
        if State == Par.State_Undefined,
            error('Undefined track');
        elseif any(State == (Par.State_FisrtInit:Par.State_LastInit-1)),    % initialization states
            NextState = State + 1;                          % next initialization state
        elseif State == Par.State_LastInit,                 % last initialization state
            NextState = Par.State_Track;                    % next initialization state
        elseif any(State == (Par.State_FirstCoast:Par.State_LastCoast)),    % coast mode states
            NextState = Par.State_Track;                    % return to track state
        else                                                % track state
            NextState = Par.State_Track;                    % stay in this state
        end;
        
        % track velocity initalization
        if (State == Par.State_FisrtInit) & 1,
            
            % observed data is OK !!! his is not a real PDAF update !!!
            y = mean(DataList(:,DataInd),2);
 
           % x is initialized with current measurements
            yprev = H*x;
            veloc = (y - yprev)./Par.dT;
            % ObservInd + 1 gives us velocity index
            x(TrackList{i}.TrackObj.ObservInd+1) = veloc;
            
        end;
        
    end;

    
    %------------------------------------------------
    % make prediction
    xpred   = F*x;
    Ppred   = F*P*F' + Q;
    
    ypred   = H*xpred;
    S       = H*Ppred*H' + R;
    Sinv    = inv(S);

    %------------------------------------------------
    % check the data
    if isempty(DataInd),
        % use predicted data
        y = H*xpred;
    else
        % associated data
        if Par.UsePDAF, 
           y = DataList(:,DataInd);
        else
            y = mean(DataList(:,DataInd),2); 
        end;
    end;
    [ObsDim,ValidDataPointNum] = size(y);
    
    %----------------------------------
    % basic kalman filter
    e       = y - ypred(:,ones(1,ValidDataPointNum)); % error (innovation) for each sample
    W       = Ppred*H'*Sinv;                    % Kalman gain matrix
    Pc      = (eye(size(F,1)) - W*H)*Ppred;

    
    %------------------------------------------------
    % compute association probabilities
    if Par.UsePDAF,
        
        loglik  = sum((e'*Sinv).*e',2); % volume computation for the vector     
        betta   = zeros(1,ValidDataPointNum+1); % the last one is zero
        % The next is OK
        betta(1:ValidDataPointNum) = exp(-.5*loglik); % a
        betta(ValidDataPointNum+1) = (1-PG*PD)/PD*2*ValidDataPointNum/Par.GateLevel*sqrt(det(S)); % by Tracking and Data Association
        %betta(vlen+1) = vlen*(1-PG*PD)/PD/(pi*gamma*sqrt(det(S))); % by Multitarget-Multisensor Tracking
        betta   = betta./sum(betta);
	
        %------------------------------------------------
        % update
        etot    = e*betta(1:ValidDataPointNum)';
        Pgag    = W*((e.*betta(ones(ObsDim,1),1:ValidDataPointNum))*e' - etot*etot')*W';
        Pnew    = betta(ValidDataPointNum+1)*Ppred + (1-betta(ValidDataPointNum+1))*Pc + Pgag;
    else
        etot    = e;
        Pnew    = Pc;
    end;
        
    xnew    = xpred + W*etot;

    
   % This is not a real PDAF update
    TrackList{i}.TrackObj.x     = xnew;
    TrackList{i}.TrackObj.P     = Pnew;
    TrackList{i}.TrackObj.State = NextState;
    
    % history update
    Hist(:,2:end)               = Hist(:,1:end-1);
    Hist(:,1)                   = x;
    TrackList{i}.TrackObj.Hist  = Hist;
    
    % life time update
    TrackList{i}.TrackObj.LifeTime = TrackList{i}.TrackObj.LifeTime + 1;
    
    % likelihood update
    %if NextState == Par.State_Track, % data is found and OK
    if any(NextState == (Par.State_FirstCoast:Par.State_LastCoast)), % no data is found - make worth update
        TrackList{i}.TrackObj.LogLike  = (1-alpha)*TrackList{i}.TrackObj.LogLike + alpha*Par.GateLevel;
    else                             % init state or coast log likelihood is increased
        TrackList{i}.TrackObj.LogLike  = (1-alpha)*TrackList{i}.TrackObj.LogLike + alpha*etot'*Sinv*etot;
        %TrackList{i}.TrackObj.LogLike  = (1+alpha)*TrackList{i}.TrackObj.LogLike;
    end;
end;    % track loop
