function [TrackList] = Track_InitConfDel(TrackList,DataList,ValidationMatrix,bettaNTFA, betta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track_Maintenance - Performs track initiation, confirmation and deletion 
% Input:
%   TrackList        - List of Tracks
%   ValidationMatrix - Matrix showing all valid data associations
% Output:
%   TrackList    - Updated list of Tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TrackNum    = size(TrackList,2);
    PointNum    = size(DataList,2);
    ObsDim      = size(DataList,1);
    P_TM = 0.1;                         % Track miss probability (from Blackman & Popoli)
    P_FC = 0.00001;                     % False confirm probability (from Blackman & Popoli)
    P_D = 0.8;                          % Probability of detection
    P_G = 0.9;                          % Probability of gating
    THD = -6.9;                         % Deletion threshold diff (from Blackman & Popoli)
    Vmax = 0.4;
    
    % High and low thresholds for confirmation and deletion
    %  of tentative tracks
    gamma_low = log(P_TM/(1-P_FC));
    gamma_high = log((1-P_TM)/P_FC);
    
    % Kalman Parameters
    n=4;      %number of state
    q=0.01;    %std of process 
    r=0.25;    %std of measurement
    s.Q=[1^3/3, 0, 1^2/2, 0;  0, 1^3/3, 0, 1^2/2; 1^2/2, 0, 1, 0; 0, 1^2/2, 0, 1]*10*q^2; % covariance of process
    s.R=r^2*eye(n/2);        % covariance of measurement  
    s.sys=(@(x)[x(1)+ x(3); x(2)+x(4); x(3); x(4)]);  % assuming measurements arrive 1 per sec
    s.obs=@(x)[x(1);x(2)];                               % measurement equation                                % initial state
    s.x_init = [];
    s.P_init = [];

    
    % Two-point difference initiation
    %for i = 1:TrackNum
    %    s.x_init(:,i)=[DataList(1,i,2); DataList(2,i,2); DataList(1,i,2)-DataList(1,i,1); DataList(2,i,2)-DataList(2,i,1)]; %initial state
    %end
    %s.P_init = [q^2, 0, q^2, 0;
    %            0, q^2, 0, q^2;
    %            q^2, 0, 2*q^2, 0;
    %            0, q^2, 0, 2*q^2];                               % initial state covraiance

    % Single-point initiatio


    % Process and Observation handles and covariances
    TrackObj.sys          = s.sys;
    TrackObj.obs          = s.obs;
    TrackObj.Q          = s.Q;
    TrackObj.R          = s.R;

    % initial covariMeasInd, TrackIndance assumed same for all tracks
    TrackObj.P          = s.P_init;
    
    % if the ValidationMatrix is empty (signals overall initiation step)
    if numel(ValidationMatrix)==0
        UnassocMeasInd = [1:PointNum];
    else
        % Get indices of all rows (measurements) where there exist no possible association
        UnassocMeasInd = find(all(ValidationMatrix==0,2));
    end

    % Update LPR and status of existing existing tracks
    %ConfirmedTracks = find(cell2mat(cellfun(@(sas)sas.Status, TrackList.TrackObj, 'uni', false))=='Confirmed');
    DeletedTracks = [];
    for i=1:TrackNum
        TrackInd = i;
        PossibleAssocMeas = find(ValidationMatrix(:,TrackInd));

        % Update LPR
        if isempty(PossibleAssocMeas)
            TrackList{TrackInd}.TrackObj.LPR = TrackList{TrackInd}.TrackObj.LPR + log(1-P_D*P_G);
        else
            for j=1:numel(PossibleAssocMeas)
                TrackList{TrackInd}.TrackObj.LPR = TrackList{TrackInd}.TrackObj.LPR + betta(TrackInd,PossibleAssocMeas(j))*log(P_D*mvnpdf(DataList(:,PossibleAssocMeas(j)),TrackList{TrackInd}.TrackObj.z_pred, TrackList{TrackInd}.TrackObj.S)/bettaNTFA);
            end

            % Update maximum LPR
            if (TrackList{TrackInd}.TrackObj.LPR>TrackList{TrackInd}.TrackObj.LPR_max)
                TrackList{TrackInd}.TrackObj.LPR_max = TrackList{TrackInd}.TrackObj.LPR;
            end
        end

        % Check against thresholds
        if (TrackList{TrackInd}.TrackObj.Status==0) % Tentative tracks
            if TrackList{TrackInd}.TrackObj.LPR>gamma_high
                TrackList{TrackInd}.TrackObj.Status = 1; % Confirm track
            elseif (TrackList{TrackInd}.TrackObj.LPR<gamma_low)
                DeletedTracks(end+1) = TrackInd;
                %TrackList{TrackInd} = []; % Delete Track
            end
        else % Confirmed tracks
            if (TrackList{TrackInd}.TrackObj.LPR < TrackList{TrackInd}.TrackObj.LPR_max+THD)
                DeletedTracks(end+1) = TrackInd;
                %TrackList{TrackInd} = []; % Delete Track
            end
        end  
    end
    
    for i=1:numel(DeletedTracks)
        TrackList(DeletedTracks(i)) = [];
    end
    
%     if (~isempty(TrackList))
%         % Remove cells of deleted tracks from cell array
%         TrackList(~cellfun('isempty',TrackList));
%     end
        
    % Initiate new tracks
    for i=1:numel(UnassocMeasInd)
        MeasInd = UnassocMeasInd(i);
        s.x = [DataList(1,MeasInd); DataList(2,MeasInd); 0; 0]; %initial state
        s.P = diag([q^2, q^2, (Vmax^2/3), (Vmax^2/3)]);
        s.Status = 0;
        s.LPR = 1;
        s.LPR_max = s.LPR;
        TrackList{end+1}.TrackObj = s; 
    end
end