function DataList = gen_obs_cluttered_multi(TrackNum, x_true, y_true)
    r = 0.5
    r_clutter = 2;
    DataList = [];
    for i=1:size(x_true,1)

        for j=1:TrackNum
            if (x_true(i,j)~=0 && y_true(i,j)~=0)
                DataList(:,j,i) = [ x_true(i,j)+normrnd(0,r^2), y_true(i,j)+normrnd(0,r^2)];
            end
        end
        
        % FA for target 1 (always present)
        for j=TrackNum+1:10
             DataList(:,j,i) = [ x_true(i,1)+normrnd(0,r_clutter^2), y_true(i,1)+normrnd(0,r_clutter^2)];
        end
        
        % FA for target 2 
        if (TrackNum>=2)
            for j=11:17
                 DataList(:,j,i) = [ x_true(i,2)+normrnd(0,r_clutter^2), y_true(i,2)+normrnd(0,r_clutter^2)];
            end
        end
        
        % FA for target 3
        if (TrackNum>=3)
            for j=18:25
                 DataList(:,j,i) = [ x_true(i,3)+normrnd(0,r_clutter^2), y_true(i,3)+normrnd(0,r_clutter^2)];
            end
        end
    end
end