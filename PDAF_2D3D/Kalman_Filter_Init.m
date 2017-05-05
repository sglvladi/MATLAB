function [F,H,Q,R,ObservInd,initx,initP] = Kalman_Filter_Init(dT,ModelDim,ProbDim,StateVar,ObserVar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kalman_Filter_Init - initializes all the Kalman filter
% matrix with respect to dimensions
% Input:
%   dT          - time sampling interval
%   ProbDim     - dimensionality of the problem (x, x-y or x-y-z)
%   ModelDim    - constant velocity or constant acceleration
%   StateVar    - maneuver parameter - sqrt(StateVar)
%   ObserVar    - noise variance - sqrt(ObserVar)
% Output:
%   F,H,Q,R     - kalman structure
%   ObservInd   - index of observable variables
%   initx,initP - initial x and P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters
if nargin < 5, ObserVar = 0.01; end;
if nargin < 4, StateVar = 1; end;% maneuver parameter - sqrt(StateVar)
if nargin < 3, ProbDim = 2; end;
if nargin < 2, ModelDim = 2; end;
if nargin < 1, error('Non sufficient input parameters'); end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model intialization

%StateVar = dT * 1; % maneuver parameter - sqrt(StateVar)
%ObserVar = dT * 0.01; %.005;

% possible model
% A = [0 1 0;0 0 1;0 0 -.1]; F = eye(3)+dt*A+.5*(dt*A)^2+1/6*(dt*A)^3
switch ModelDim,
case 2,
    Ftmp = [1 dT;0 1];
    Htmp = [1 0];
    Qtmp = [dT^4/4 dT^3/2; dT^3/2 dT^2]*StateVar; % Multisensor track. of maneuvr. target in clutter
    %Qtmp = [(dT^2)/3 dT/2; dT/2 1]*dT*StateVar;
    Rtmp = 1*ObserVar;
    xtmp = [0;0];
    Ptmp = Qtmp*10;
    
case 3,
    Ftmp = [1 dT dT^2/2;0 1 dT;0 0 1];
    Htmp = [1 0 0];
    Qtmp = [(dT^4)/20 (dT^3)/8 (dT^2)/6;... % from Kalman_tutorial_paper
            (dT^3)/8  (dT^2)/3  dT/2   ;...
            (dT^2)/6   dT/2       1    ]*dT*StateVar;
%     Qtmp = [(dT^4)/4 (dT^3)/2 (dT^2)/2;... % Multisensor track. of maneuvr. target in clutter
%             (dT^3)/2  dT^2      dT    ;...
%             (dT^2)/2   dT        1    ]*StateVar;
    Rtmp = 1*ObserVar;
    xtmp = [0;0;0];
    Ptmp = Qtmp*100;
    
otherwise
    error('Unsupported model dimensions.') 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem dimensionality
switch ProbDim,
case 2,             % 2D stacking
    F       = [Ftmp zeros(size(Ftmp)); zeros(size(Ftmp)) Ftmp];
    H       = [Htmp zeros(size(Htmp)); zeros(size(Htmp)) Htmp]; 
    Q       = [Qtmp zeros(size(Qtmp)); zeros(size(Qtmp)) Qtmp]; 
    R       = [Rtmp zeros(size(Rtmp)); zeros(size(Rtmp)) Rtmp]; 
    initP   = [Ptmp zeros(size(Ptmp)); zeros(size(Ptmp)) Ptmp]; 
    initx   = [xtmp; xtmp]; 
    ObservInd = [1 size(Ftmp,1)+1];
case 3,             % 3D stacking
    F       = kron(eye(3),Ftmp);
    H       = kron(eye(3),Htmp); 
    Q       = kron(eye(3),Qtmp); 
    R       = kron(eye(3),Rtmp);   
    initP   = kron(eye(3),Ptmp); 
    initx   = kron(ones(3,1),xtmp);
    ObservInd = [1 size(Ftmp,1)+1 2*size(Ftmp,1)];    
otherwise
    error('Unsupported problem dimensions.')    
end;