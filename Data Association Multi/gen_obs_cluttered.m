r = 0.5
r_clutter = 2;
DataList = [];
for i=1:size(x_true)
    DataList(:,1,i) = [ x_true(i)+r*normrnd(0,r^2), y_true(i)+r*normrnd(0,r^2)];
    for j=2:6
        DataList(:,j,i) = [ x_true(i)+r_clutter*normrnd(0,r_clutter^2), y_true(i)+r_clutter*normrnd(0,r_clutter^2)];
    end
end