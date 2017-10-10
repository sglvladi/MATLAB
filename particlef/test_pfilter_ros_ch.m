%% clear memory, screen, and close all figures

%% Process equation x[k] = sys(k, x[k-1], u[k]);
nx = 5;  % number of states
sys = @(k, xkm1, uk) [xkm1(1)+k*xkm1(3)*cos(xkm1(4)); xkm1(2)+k*xkm1(3)*sin(xkm1(4)); xkm1(3) + uk(3); xkm1(4) + uk(4)]; % (returns column vector)

%% Observation equation y[k] = obs(k, x[k], v[k]);
ny = 2;                                           % number of observations
obs = @(k, xk, vk) [xk(1)+vk(1); xk(2)+vk(2)];                  % (returns column vector)

%% PDF of process noise and noise generator function
nu = 5;                                           % size of the vector of process noise
sigma_u = 0.01;
cov_u = [1^3/3, 0, 1^2/2, 0, 0;  0, 1^3/3, 0, 1^2/2, 0; 1^2/2, 0, 1, 0, 0; 0, 1^2/2, 0, 1, 0; 0 0 0 0 0]*10*sigma_u^2;
p_sys_noise   = @(u) mvnpdf(u, zeros(1, nu), cov_u);
gen_sys_noise = @(u) mvnrnd(zeros(1, nu), cov_u);         % sample from p_sys_noise (returns column vector)

%% PDF of observation noise and noise generator function
nv = 2;                                           % size of the vector of observation noise
sigma_v = 0.5;
cov_v = sigma_v^2*eye(nv);
p_obs_noise   = @(v) mvnpdf(v, zeros(1, nv), cov_v);
gen_obs_noise = @(v) mvnrnd(zeros(1, nv), cov_v);         % sample from p_obs_noise (returns column vector)

%% Initial PDF
% p_x0 = @(x) normpdf(x, 0,sqrt(10));             % initial pdf
gen_x0 = @(x) mvnrnd([obs_x(2); obs_y(2); obs_x(2)-obs_x(1); obs_y(2)-obs_y(1);0],cov_u);               % sample from p_x0 (returns column vector)

%% Transition prior PDF p(x[k] | x[k-1])
% (under the suposition of additive process noise)
% p_xk_given_xkm1 = @(k, xk, xkm1) p_sys_noise(xk - sys(k, xkm1, 0));

%% 
% (under the suposition of additive process noise)
p_yk_given_xk = @(k, yk, xk) p_obs_noise((yk - obs(k, xk, zeros(1, nv)))');

%% Number of time steps
T = 340;

%% Separate memory space
%x = zeros(nx,T);  y = zeros(ny,T);
%u = zeros(nu,T);  v = zeros(nv,T);

%% Simulate system
xh0 = [obs_x(2); obs_y(2); obs_x(2)-obs_x(1); obs_y(2)-obs_y(1); 0]                                  % initial state
%u(:,1) = gen_sys_noise(sigma_u)';                               % initial process noise
%v(:,1) = gen_obs_noise(sigma_v)';          % initial observation noise
% x(:,1) = xh0;
% y(:,1) = obs(1, xh0, v(:,1));
% for k = 2:T
%    % here we are basically sampling from p_xk_given_xkm1 and from p_yk_given_xk
%    %u(:,k) = gen_sys_noise();              % simulate process noise
%    v(:,k) = gen_obs_noise();              % simulate observation noise
%    x(:,k) = sys(k, x(:,k-1), u(:,k));     % simulate state
%    y(:,k) = obs(k, x(:,k),   v(:,k));     % simulate observation
% end
%x = [x_true, y_true]';
y = [obs_x'; obs_y'];
%% Separate memory
xh = zeros(nx, T); xh(:,1) = xh0;
yh = zeros(ny, T); yh(:,1) = obs(1, xh0, zeros(1, nv));

pf.k               = 1;                   % initial iteration number
pf.Np              = 1000;                 % number of particles
%pf.w               = zeros(pf.Np, T);     % weights
pf.particles       = zeros(nx, pf.Np, T); % particles
pf.gen_x0          = gen_x0;              % function for sampling from initial pdf p_x0
pf.obs             = p_yk_given_xk;       % function of the observation likelihood PDF p(y[k] | x[k])
pf.sys_noise       = gen_sys_noise;       % function for generating system noise
%pf.p_x0 = p_x0;                          % initial prior PDF p(x[0])
%pf.p_xk_given_ xkm1 = p_xk_given_xkm1;   % transition prior PDF p(x[k] | x[k-1])
pf.xhk = xh0;
pf.sys = sys;
pf.resampling_strategy = 'systematic_resampling';
my_pf = ParticleFilter(pf);

figure
%% Estimate state
for k = 2:T
   clf;
   fprintf('Iteration = %d/%d\n',k,T);
   % state estimation
   my_pf.pf.k = Dt;
   my_pf.pf.z = y(:,k);
   %[xh(:,k), pf] = particle_filter(sys, y(:,k), pf, 'multinomial_resampling');
   %[xh(:,k), pf] = particle_filter(sys, y(:,k), pf, 'systematic_resampling');   
   my_pf.pf = my_pf.Iterate(my_pf.pf);
   xh(:,k) = my_pf.pf.xhk(:,k);
   % filtered observation
   yh(:,k) = obs(k, xh(:,k), zeros(1, nv));
   pV_err(1,k) = sqrt((xh(1,k) - x(1,k))^2 + (xh(2,k)-x(2,k))^2);
%    hold on
%    sort_w = sort(my_pf.pf.w(:,k))
%    for i=1:pf.Np
%        if my_pf.pf.w(i,k)<=sort_w(10)
%            plot(my_pf.pf.particles(1,i,k), my_pf.pf.particles(2,i,k),'r.')
%        elseif my_pf.pf.w(i,k)<=sort_w(500)
%            plot(my_pf.pf.particles(1,i,k), my_pf.pf.particles(2,i,k),'y.')
%        else
%            plot(my_pf.pf.particles(1,i,k), my_pf.pf.particles(2,i,k),'c.')
%        end
%    end
   %plot(my_pf.pf.particles(1,:,k), my_pf.pf.particles(2,:,k),'c.')
   hold on;
   h2 = plot(x(1,1:k),x(2,1:k),'b','LineWidth',1);
   h3 = plot(xh(1,1:k),xh(2,1:k),'m','LineWidth',1);
   h4 = plot(y(1,1:k), y(2,1:k),'g.','LineWidth',1);
   h5 = plot(xh(1,1:k),xh(2,1:k),'m.','LineWidth',3);
   axis([0 10 0 10])
   legend([h2 h3 h4],'state','mean of estimated state','measurements');
   title('State vs estimated state by the particle filter vs particle paths','FontSize',14);
    pause(0.1)

end

%% Compute RMSE
err = (xh(1:2,:) - x(:,1:T).*(xh(1:2,:) - x(:,1:T)));
RMSE = sqrt(sum(err,2)/T)

%% Make plots of the evolution of the density
figure
hold on;
xi = 1:T;
yi = -25:0.25:25;
[xx,yy] = meshgrid(xi,yi);
den = zeros(size(xx));
xhmode = zeros(size(xh));

for i = xi
   % for each time step perform a kernel density estimation
   den(:,i) = ksdensity(my_pf.pf.particles(1,:,i), yi,'kernel','epanechnikov');
   [~, idx] = max(den(:,i));

   % estimate the mode of the density
   xhmode(i) = yi(idx);
   plot3(repmat(xi(i),length(yi),1), yi', den(:,i));
end
view(3);
box on;
title('Evolution of the state density','FontSize',14)

figure
mesh(xx,yy,den);   
title('Evolution of the state density','FontSize',14)

%% plot of the state vs estimated state by the particle filter vs particle paths
figure
hold on;
%h1 = plot(1:T,squeeze(my_pf.pf.particles(2,:,:)),'y');
h2 = plot(x(1,:),x(2,:),'b','LineWidth',1);
h3 = plot(xh(1,:),xh(2,:),'r','LineWidth',1);
h4 = plot(y(1,:), y(2,:),'g.','LineWidth',1);
legend([h2 h3 h4], 'state','mean of estimated state','measurements');
title('State vs estimated state by the particle filter vs particle paths','FontSize',14);

%% plot of the observation vs filtered observation by the particle filter
figure
plot(1:T,y(1,:),'b', 1:T,yh(1,:),'r');
legend('observation','filtered observation');
title('Observation vs filtered observation by the particle filter','FontSize',14);

% %% Plot results
% figure
% for k=1:4                                 % plot results
%     subplot(4,1,k)
% %     figure
% %     hold on
%     plot( 1:T, squeeze(my_pf.pf.particles(k,:,:)), 'y', 1:T, x(k,:), 'k--', 1:T, xh(k,:), 'b', 1:T, y(k,:), 'r.')
% %     for i = 1:N
% %         hold on
% %         error_ellipse('C', blkdiag(PV_ukf{i}(1,1),1), 'mu', [i, xV_ukf(k,i)], 'style', 'r--')
% %         hold on
% %         error_ellipse('C', blkdiag(PV_ekf{i}(1,1),1), 'mu', [i, xV_ekf(k,i)], 'style', '--')
% %     end
%     str = sprintf('Patricle Filter Estimated State X(%d)',k);
%     title(str)
%     legend('Real', 'Particle paths', 'Estimated Mean', 'Meas');
% end
% return;