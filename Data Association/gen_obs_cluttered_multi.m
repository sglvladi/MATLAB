function [DataList, x, y] = gen_obs_cluttered_multi(TrackNum, x_true, y_true, q , q_clutter, No_of_clutter)
%     q = 0.5
%     q_clutter = 8;
    DataList = [];
    k=1;
    x = [];
    y = [];
    for i=1:2:size(x_true,1)
        for j=1:TrackNum
            x(k,j) = x_true(i,j);
            y(k,j) = y_true(i,j);
            DataList(:,j,k) = [ x_true(i,j)+normrnd(0,q), y_true(i,j)+normrnd(0,q)];
        end
        
        for j=TrackNum+1:No_of_clutter
            DataList(:,j,k) = [ x_true(i,1)+normrnd(0,q_clutter), y_true(i,1)+normrnd(0,q_clutter)];
        end
        k=k+1;
    end
end