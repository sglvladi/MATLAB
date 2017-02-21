function TrackShow = Structure_PDAF_Show(TrackShow,Par,TrackList,DataList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure_PDAF_Show - draws the data and tracks for each time interval
% Input:
%   TrackList    - kalman structure list and more
%   DataList     - 2 x DataPointsNum contains the relevant data for time t
%   Par          - parameters
% Output:
%   TrackShow    - structure to hold handles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First time init
if nargin < 3 || isempty(TrackShow),
    % do init
    TrackShow  = ShowInit(Par);
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%persistent ShowFigNum;
ShowFigNum          = TrackShow.ShowFigNum;            % show figure number
%AxisSc      = [Par.Y1Bounds Par.Y1Bounds]*1.1+[-.1 0 -.1 0];
SmallShift          = TrackShow.SmallShift;     % used for text display
%NumSigma            = TrackShow.NumSigma;%sqrt(Par.GateLevel); % number of sigmas


% handles
DataPlotHand        =  TrackShow.DataPlotHand;
TrackPlotHand       = TrackShow.TrackPlotHand;
TrackNumHand        = TrackShow.TrackNumHand;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show points
%DLen = size(DataList,2);
figure(ShowFigNum), hold on;

if ~isempty(DataPlotHand),
    delete(DataPlotHand);
end;
switch Par.ProblemDim,
     case 2,
        DataPlotHand        = plot(DataList(1,:), DataList(2,:),'b.'); 
     case 3,
        DataPlotHand        = plot3(DataList(1,:), DataList(2,:),DataList(3,:),'b.'); 
 end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show tracking objects   
% loop throught all tracks
TrackNum            = length(TrackList);
ValidStatesForShow  = [Par.State_Track Par.State_FirstCoast:Par.State_LastCoast];

    
 

% print state development
TrackNumStr = [];
TrackStateStr = [];
% 
% if isempty(TrackPlotHand),
%     
%      for i=1:TrackNum,
%          y = TrackList{i}.TrackObj.H * TrackList{i}.TrackObj.x;
%          H = TrackList{i}.TrackObj.H;
%          P = TrackList{i}.TrackObj.P;
%          R = TrackList{i}.TrackObj.R;
%          
%          S = H*P*H' + R; 
%          [u,sing,v] = svd(S);
%          %elipse = u*sing*circle;
%          elipse = u*sqrt(sing)*circle;
%          TrackPlotHand(i) = plot(elipse(1,:)+y(1), elipse(2,:)+y(2),'r'); 
%          TrackNumHand(i) = text(y(1)+SmallShift,y(2),num2str(i),'FontSize',8);
%          
%          TrackNumStr = [TrackNumStr sprintf(' %2.0d',i)];
%          
%      end;
% 
%     % one time print
%     if ~Par.ShowOn, fprintf('%s\n',TrackNumStr); end;
%    
%  end;
 
 if ~isempty(TrackPlotHand)
 delete(TrackPlotHand);
 delete(TrackNumHand);
 end
  
 for i=1:TrackNum,
     
     
     y = TrackList{i}.TrackObj.H * TrackList{i}.TrackObj.x;
     H = TrackList{i}.TrackObj.H;
     P = TrackList{i}.TrackObj.P;
     R = TrackList{i}.TrackObj.R;
     
     S = H*P*H' + R;          
     %S

     [u,sing,v] = svd(S);
     elipse     = u*sqrt(sing)*TrackShow.Ball;
      %elipse = u*sing*circle;
     
     if ~any(TrackList{i}.TrackObj.State == ValidStatesForShow),
        y = nan(Par.ProblemDim,1); % to prevent show
     end;

     switch Par.ProblemDim,
         case 2,
             TrackPlotHand(i) = plot(elipse(1,:)+y(1), elipse(2,:)+y(2),'r'); 
             TrackNumHand(i) = text(y(1)+SmallShift,y(2),num2str(i),'FontSize',8);
         case 3,
             TrackPlotHand(i) = plot3(elipse(1,:)+y(1), elipse(2,:)+y(2),elipse(3,:)+y(3),'r'); 
             TrackNumHand(i) = text(y(1)+SmallShift,y(2),y(3)+SmallShift,num2str(i),'FontSize',8);
     end
     
     TrackStateStr = [TrackStateStr sprintf(' %2.0d',TrackList{i}.TrackObj.State)];
     %TrackStateStr = [TrackStateStr sprintf(' %2.0i',round(TrackList{i}.TrackObj.LogLike*100))];
     
 end;

 if ~Par.ShowOn, fprintf('%s\n',TrackStateStr); end;
       
drawnow;
hold off;

% save
TrackShow.DataPlotHand        = DataPlotHand;
TrackShow.TrackPlotHand       = TrackPlotHand;
TrackShow.TrackNumHand        = TrackNumHand;


return




%%%%%%%%
%%%%%%
% Internal Function
 function TrackShow  = ShowInit(Par)
 % Init params for show
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%persistent ShowFigNum;
TrackShow.ShowFigNum  = 15;            % show figure number
TrackShow.SmallShift  = 5e-3;     % used for text display



NumSigma    = sqrt(Par.GateLevel); % number of sigmas

figure(TrackShow.ShowFigNum), clf, hold on;
switch Par.ProblemDim,
    case 2,
        TrackShow.AxisSc        = [Par.Y1Bounds Par.Y2Bounds]*1.1+[-.1 0 -.1 0];
        TrackShow.DataPlotHand  = [];
        TrackShow.TrackPlotHand = [];
        TrackShow.TrackNumHand  = [];
        
        TrackShow.DataPlotHand   = plot(0, 0, 'b.'); 
        axis(TrackShow.AxisSc);
        title('Data Points & Tracks')
        xlabel('X1'),ylabel('X2')
        
        % for graphics
        t                       = 0:.1:2*pi+.1;
        TrackShow.Ball          = [cos(t);sin(t)]*NumSigma; %
       
        
    case 3,
        TrackShow.AxisSc        = [Par.Y1Bounds Par.Y2Bounds Par.Y3Bounds]*1.1+[-.1 0 -.1 0 -.1 0];
        TrackShow.DataPlotHand  = [];
        TrackShow.TrackPlotHand = [];
        TrackShow.TrackNumHand  = [];
        TrackShow.DataPlotHand   = plot3(0, 0, 0,'b.'); 
        axis(TrackShow.AxisSc);
        title('Data Points & Tracks')
        xlabel('X1'),ylabel('X2'),zlabel('X3')
        view(17,34)
  
        % for graphics
        [X,Y,Z]                 = sphere(10);
        TrackShow.Ball          = [X(:)';Y(:)';Z(:)']*NumSigma; %

        
    otherwise
        error('Par.ProblemDim')
end

box on;     
hold off;
     
 return

