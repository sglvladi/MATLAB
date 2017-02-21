% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Test EKF vs UKF example                                      %
% ------------------------------------------------------------ %
% Author: Lyudmil Vladimirov                                   %
% University of Livepool                                       %
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
clearvars ukf_x ukf_P;

%% Get size of measurement table
[m, n] = size(meas);

%% Initialise variables:
clear s
s.x = zeros(4,1);       % state vector
f_pre=@(Dt)(@(x)[x(1)+Dt*x(3);x(2)+Dt*x(4);x(3);x(4)]); % Time-dependant Process function
s.sys = f_pre(Dt);        % Process function
s.Q = zeros(4,4);       % Process noise
sigma = 0;              % sigma
s.obs = @(x)[x(1);x(2)];  % Measurement function (Likelihood)
s.R = [sigma^2 0; 0 sigma^2];       % Measurement noise
s.B = 0;                % No control inputs               
s.u = 0;                %       >>

%% Set initial time
t_s = 2;

%% Populate initial state vector by 2 point initialization
Dt = second(Timestamp(t_s,1)) -  second(Timestamp(t_s-1,1));
if Dt~=0
    s.x = [ meas(2,1); meas(2,2); (meas(2,1)-meas(1,1))*Dt; (meas(2,2)-meas(1,2))*Dt]; 
else
    while Dt==0
        t_s = t_s + 1;
        Dt = second(Timestamp(t_s,1)) -  second(Timestamp(t_s-1,1));
        s.x = [ meas(t_s,1); meas(t_s,2); (meas(t_s,1)-meas(t_s-1,1))*Dt; (meas(t_s,2)-meas(t_s-1,2))*Dt]; 
    end
end

if sigma~=0
    s.P = [sigma^2, 0, sigma^2/Dt, 0; 0, sigma^2, 0, sigma^2/Dt; sigma^2/2, 0, 2*sigma^2/Dt^2, 0; 0, sigma^2/2, 0, 2*sigma^2/Dt^2];
else
    s.P = eye(4);
end

%% Initialise EKF and UKF
kf = EKalmanFilter(s);
ukf = UKalmanFilter(s, 1, -1, 2);

%% Dummy variables used to store data for display
x = [];
P = [];
ukf_x = [];
ukf_P = [];
z = [];
smooth = kf.s;
x_ukf = kf.s(end).x;
P_ukf = kf.s(end).P;

xk_ukf = x_ukf;
Pk_ukf = P_ukf;
%% Main loop
for t=t_s+1:2262
    
    % Compute Dt
    datestr(Timestamp(t,1))
    t1 = datevec(datestr(Timestamp(t,1)),'dd-mm-yyyy HH:MM:SS');
    t2 = datevec(datestr(Timestamp(t-1,1)),'dd-mm-yyyy HH:MM:SS');
    Dt = etime(t1,t2);
    
    % Only iterate through filter if Dt different than 0
    if Dt~=0
        xkm1_ukf = xk_ukf;
        Pkm1_ukf = Pk_ukf;
        kf.s(end).z = transpose(meas(t,:)); % Fetch next measurement
        kf.s(end).sys = f_pre(Dt);            % Compute process function given Dt
        kf.s(end).Q = [Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt];    % Compute process noise given Dt
        ukf.s(end).z = transpose(meas(t,:)); % Fetch next measurement
        ukf.s(end).sys = f_pre(Dt);            % Compute process function given Dt
        ukf.s(end).Q = [Dt^3/3, 0, Dt^2/2, 0;  0, Dt^3/3, 0, Dt^2/2; Dt^2/2, 0, Dt, 0; 0, Dt^2/2, 0, Dt];   % Compute process noise given Dt

        kf.s = kf.Iterate(kf.s);        % Run EKF 
        ukf.s = ukf.Iterate(ukf.s);     % Run UKF
        %sys_pre = @(Dt)(@(x,y)[x(1)+Dt*x(3);x(2)+Dt*x(4);x(3);x(4)]);
        %obs = @(x,y)[x(1);x(2)];  % Measurement function (Likelihood)
        %[xk_ukf, Pk_ukf] = ukf_fl(sys_pre(Dt), obs, kf.s(end).z, xkm1_ukf, Pkm1_ukf, kf.s(end).Q, kf.s(end).R);
        
    end
    
    % Store data in dummy variables for display
    smooth(end+1) = kf.s;
    x =[x , kf.s.x];
    P =[P , kf.s.P(1,1)];
    ukf_x = [ukf_x, ukf.s.x]; %xk_ukf
    ukf_P = [ukf_P, ukf.s.P(1,1)]; %Pk_ukf(1,1)
    z = [z, kf.s.z];   
end

%% Plot measurement data:
x=x';
P=P';
ukf_x = ukf_x';
ukf_P = ukf_P';
z = z';
kf = kf.Smooth(smooth);
figure
hold on
grid on
subplot(2,1,1)
hold on
hz=plot(z(:,1),z(:,2),'r.');
% plot a-posteriori state estimates:
%x = transpose([kf.s(2:end).x]);
x_smooth = kf.x_smooth';
hk=plot(x(:,1),x(:,2),'bo-');
hk=plot(ukf_x(:,1),ukf_x(:,2),'go-');
title('EKF vs UKF estimated state sequence')
legend('GPS', 'EKF', 'UKF');

subplot(2,1,2)
hold on
hk = plot(P(:,1),'bo-');
hk = plot(ukf_P(:,1),'go-');
%x = xV';
%hy=plot(x_smooth(:,1), x_smooth(:,2), 'g-');
%hy=quiver3(s.x(:,1),s.x(:,2),s.x(:,3),s.x(:,4));
%ht=plot(tru,'g-');
%legend([hz hk ht],'observations','Kalman output','true voltage',0)
title('EKF vs UKF estimated covariance')
legend('EKF', 'UKF');
