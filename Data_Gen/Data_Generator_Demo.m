dat_gen = Data_Generator();
dat_gen = dat_gen.setProperties(4, 2, @(x,Dt) [x(1)+x(3)*0.1*Dt; x(2)+x(4)*0.1*Dt; x(3); x(4)], 5, @(x)[x(1); x(2)], 1, 10, 200);
[data] = dat_gen.Generate();
figure
for i = 1:10
    plot(data{i}.trackObj.x_true(1,:), data{1}.trackObj.x_true(2,:));
    hold on
end