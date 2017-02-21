function [DataList, x, y] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, q , q_clutter, lambda)
%     q = 0.5
%     q_clutter = 8;
    DataList = {};
    k=1;
    x = [];
    y = [];
    for i=1:1:size(x_true,1)
        for j=1:TrackNum
            x(k,j) = x_true(i,j);
            y(k,j) = y_true(i,j);
            DataList{k}(:,j) = [ x_true(i,j)+normrnd(0,q), y_true(i,j)+normrnd(0,q)];
        end
        No_of_clutter = poissrnd(lambda);
        for j=TrackNum+1:No_of_clutter
            DataList{k}(:,j) = [unifrnd(0,25), unifrnd(0,20)];
        end
        k=k+1;
    end
end