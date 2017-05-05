function TrackList = PDAF_Update(TrackList,DataList,Validation_matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Track_Update - performs track update
% Input:
%   TrackList    - List of all target tracks at time k(TrackObj's)
%   DataList     - List of all measurements at time k
%   EKF          - Kalman Filter Structure
% Output:
%   TrackList    - Updated TrackList
%
% [1] Y. Bar-Shalom, F. Daum, and J. Huang, 
%     "Probabilistic Data Association Filter", 
%     IEEE Control Systems Magazine, December 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
TrackNum    = size(TrackList,1);
% alpha       = 0.3;      % log likelihood forget factor
PG          = 0.998;      % probability of Gating
PD          = 0.8;      % probability of Detection
GateLevel   = 9;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop throught all tracks
for i=1:TrackNum,
    
    %------------------------------------------------
    % extract prediction
    x_pred       = TrackList{i}.TrackObj.x_pred;
    P_pred       = TrackList{i}.TrackObj.P_pred;
    z_pred       = TrackList{i}.TrackObj.z_pred;
    S            = TrackList{i}.TrackObj.S;
    Sinv         = inv(S);
    F            = TrackList{i}.TrackObj.F;
    H            = TrackList{i}.TrackObj.H;
    DataInd      = find(Validation_matrix(i,:));

    %------------------------------------------------
    % extract measurements
    z = DataList(:,DataInd);
    [ObsDim,ValidDataPointNum] = size(z);
    
    %----------------------------------
    % basic kalman filter
    innov_err   = z - z_pred(:,ones(1,ValidDataPointNum)); % error (innovation) for each sample
    W           = P_pred*H'*Sinv;                    % Kalman gain matrix
    Pc          = (eye(size(F,1)) - W*H)*P_pred;
    %Pc          = P_pred - W*S*W';

    
    %------------------------------------------------
    % compute association probabilities (as described in [1] for a non-parametric PDA)
    
    C   = pi; % volume of the 2-dimensional unit hypersphere     
    V_k = C*GateLevel^(ObsDim/2)*det(S)^(1/2);   % volume of the validation region 
    
    % separate memory for likelihood ratios and association prababilities
    Li_k    = zeros(1, ValidDataPointNum+1);   
    betta   = zeros(1,ValidDataPointNum+1); 
    
    % Compute likelihood ratios
    
    for j=1:ValidDataPointNum
        try
            Li_k(j) =  mvnpdf(z(:,j), z_pred, S)*PD/(ValidDataPointNum/V_k); % (ValidDataPointNum) Likelihood ratio of measurement i
        catch
            fprintf(2,'Symmetry of S has been adjusted!');
            S(2,1) = S(1,2);
            Li_k(j) =  mvnpdf(z(:,j), z_pred, S)*PD/(ValidDataPointNum/V_k); % (ValidDataPointNum) Likelihood ratio of measurement i
        end
    end
    % Compute association probabilities
    betta(1:ValidDataPointNum) = Li_k(1:ValidDataPointNum)./(1-PG*PD+sum(Li_k));
    betta(ValidDataPointNum+1) = (1-PG*PD)/(1-PG*PD+sum(Li_k));
    
    % Back-up code
    % ========================================================================>
    % loglik  = sum((innov_err'*Sinv).*innov_err',2); % volume computation for the vector
    % % The next is OK
    % betta(1:ValidDataPointNum) = exp(-.5*loglik); %a
    % betta(ValidDataPointNum+1) = (1-PG*PD)/PD*2*ValidDataPointNum/GateLevel*sqrt(det(S)); % by Tracking and Data Association
    % betta(vlen+1) = vlen*(1-PG*PD)/PD/(pi*gamma*sqrt(det(S))); % by Multitarget-Multisensor Tracking
    % betta   = betta./sum(betta);
    % ========================================================================/
    
    %------------------------------------------------
    % update
    tot_innov_err    = innov_err*betta(1:ValidDataPointNum)';
    Pgag    = W*((innov_err.*betta(ones(ObsDim,1),1:ValidDataPointNum))*innov_err' - tot_innov_err*tot_innov_err')*W';
    Pnew    = betta(ValidDataPointNum+1)*P_pred + (1-betta(ValidDataPointNum+1))*Pc + Pgag;
        
    xnew    = x_pred + W*tot_innov_err;

    
   % This is not a real PDAF update
    TrackList{i}.TrackObj.x     = xnew;
    TrackList{i}.TrackObj.P     = Pnew;
    
end;    % track loop
