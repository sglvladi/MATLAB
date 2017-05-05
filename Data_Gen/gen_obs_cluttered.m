function [x, obs_x, obs_y, obs_r, obs_phi] = gen_obs_cluttered(x_true, y_true, th_true, q, q_clutter)
    obs_x = [];
    obs_y = [];
    x = [];
    j=1;
    for i=1:2:size(x_true, 1)
        x(:,j) = [x_true(i,1); y_true(i,1); th_true(i,1)];
        noise(j,:) = mvnrnd([0,0],[q^2,0;0,q^2]);
        obs_x(j,1) = x_true(i)+noise(j,1);
        obs_y(j,1) = y_true(i)+noise(j,2);
        obs_r(j,1) = sqrt(x_true(i)^2+y_true(i)^2)+mvnrnd(0,q^2);
        obs_phi(j,1) = atan2(y_true(i),x_true(i))+mvnrnd(0,2*pi/360);
        j= j + 1;
        for k=2:6
            noise(j,:) = mvnrnd([0,0],[q_clutter^2,0;0,q_clutter^2]);
            obs_x(j,k) = x_true(i)+noise(j,1);
            obs_y(j,k) = y_true(i)+noise(j,2);
            obs_r(j,k) = sqrt(x_true(i)^2+y_true(i)^2)+mvnrnd(0,q_clutter^2);
            obs_phi(j,k) = atan2(y_true(i),x_true(i))+mvnrnd(0,2*pi/360);
        end
    end
end