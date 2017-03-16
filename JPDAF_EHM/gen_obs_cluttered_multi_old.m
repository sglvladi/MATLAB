function DataList = gen_obs_cluttered_multi_old(TrackNum, x_true, y_true)
    % Noise std
    r = 0.2;
    
    % For simplicity clutter has been normally distributed around each of
    % the targets.
    clutter_rate = 25; % Number of clutter measurements per target
    r_clutter = 2;
    
    DataList = [];
    ij = 0
    for i=1:1:size(x_true,1)
        ij=ij+1;
        for j=1:TrackNum
            if (x_true(i,j)~=0 && y_true(i,j)~=0)
                DataList(:,j,ij) = [ x_true(i,j)+normrnd(0,r^2), y_true(i,j)+normrnd(0,r^2)];
            end
        end
        
        % Clutter for target 1 (always present)
        for j=TrackNum+1:TrackNum+clutter_rate
             DataList(:,j,ij) = [ x_true(i,1)+normrnd(0,r_clutter^2), y_true(i,1)+normrnd(0,r_clutter^2)];
        end
        
        % Clutter for target 2 
        if (TrackNum>=2)
            for j=TrackNum+clutter_rate+1:TrackNum+2*clutter_rate
                 DataList(:,j,ij) = [ x_true(i,2)+normrnd(0,r_clutter^2), y_true(i,2)+normrnd(0,r_clutter^2)];
            end
        end
        
        % Clutter for target 3
        if (TrackNum>=3)
            for j=TrackNum+2*clutter_rate+1:TrackNum+3*clutter_rate
                 DataList(:,j,ij) = [ x_true(i,3)+normrnd(0,r_clutter^2), y_true(i,3)+normrnd(0,r_clutter^2)];
            end
        end
    end
end