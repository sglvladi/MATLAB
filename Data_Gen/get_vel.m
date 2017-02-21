v = [];
v(:,1) = [x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(1,1)];
for i=2:size(x_true, 1)
    v(:,i) = [x_true(i,1)-x_true(i-1,1); y_true(i,1)-y_true(i-1,1)];
end
