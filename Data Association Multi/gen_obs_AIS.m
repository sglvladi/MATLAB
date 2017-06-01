function [DataList, x, y] = gen_obs_AIS(data)
%     q = 0.5
%     q_clutter = 8;
    DataList = {};
    x = [];
    y = [];
    for i=1:size(data,1)
        DataList{i} = [data(i,1); data(i,2)];
        x(i,1) = data(i,1);
        y(i,1) = data(i,2);
    end
end