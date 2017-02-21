%function Structure_PDAF_Tracking_Demo
% Tracking a moving point in 2D plane or 3D space
% 2D : State = (x xdot y ydot). We only observe (x y).
% 3D : State = (x xdot y ydot z zdot). We only observe (x y z).

%close all,clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init parameters 
Par             = Structure_PDAF_Init_Parameters;
Par.ProblemDim  = 2;        % ProbDim     - dimensionality of the problem (x, x-y or x-y-z)
Par.TrackNum    = 3^Par.ProblemDim;        % number of trackers - space coverage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init data 
AllTheData      = Structure_PDAF_Init_Data(Par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init Kalman Filter structures
TrackList       = Structure_PDAF_Init_Track(Par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init Show structures
TrackShow       = Structure_PDAF_Show([],Par);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PDAF filtering
%FrameNum    = ;
%if Par.RecordOn, M = moviein(FrameNum); end;
for k = 1:size(AllTheData,3),
    
    % show
    if Par.ShowOn, fprintf('%3.0d : -----------------\n',k); end;
    
    % get the data for time k
    DataList    = AllTheData(:,:,k); % get the data at time 1
    
    % show
    TrackShow  = Structure_PDAF_Show(TrackShow,Par,TrackList,DataList);
    if Par.RecordOn, M(:,k) = getframe; end;


    % data-track association
    TrackList = Structure_PDAF_Association(TrackList,DataList,Par);
        
    % track update
    TrackList = Structure_PDAF_Track_Update(TrackList,DataList,Par);
    
    % track separation
    TrackList = Structure_PDAF_Track_Separation(TrackList,DataList,Par);    
    
    % start new tracks
    TrackList = Structure_PDAF_Track_Start(TrackList,DataList,Par);
    
    % recording
    %Record = Structure_PDAF_Record(Record,TrackList,DataList);

end;
if Par.RecordOn,  figure, movie(M); end;
%movie2avi(M,'Tracks6Clutter12.avi','fps',5,'quality',95)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show results
%Structure_PDAF_Show(DataList,TrackList);


