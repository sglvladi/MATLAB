q = 0.4
q_clutter = 2;
DataList = [];
for i=1:size(x_true)
    DataList(:,1,i) = [ x_true(i)+q*normrnd(0,q^2), y_true(i)+q*normrnd(0,q^2)];
    for j=2:6
        DataList(:,j,i) = [ x_true(i)+q_clutter*normrnd(0,q_clutter^2), y_true(i)+q_clutter*normrnd(0,q_clutter^2)];
    end
end