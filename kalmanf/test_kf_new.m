% SCALAR EXAMPLE (Automobile Voltimeter):

[m, n] = size(meas);
% Define the system as a constant of 12 volts:
clear s
s.x = zeros(4,1);
s.A = zeros(4,4);

% Define a process noise (stdev) of 2 volts as the car operates:
s.Q = zeros(4,4); % variance, hence stdev^2

s.P = eye(4);
% Define the voltimeter to measure the voltage itself:
s.H = [1 0 0 0; 0 1 0 0];
% Define a measurement error (stdev) of 2 volts:
s.R = [2 0; 0 2]; % variance, hence stdev^2
% Do not define any system input (control) functions:
s.B = 0;
s.u = 0;

t_s = 2;
Dt = second(Timestamp(2,1)) -  second(Timestamp(1,1));
if Dt~=0
    s.x = [ meas(2,1); meas(2,2); (meas(2,1)-meas(1,1))*Dt; (meas(2,2)-meas(1,2))*Dt]; 
end

while Dt==0
    Dt
    t_s = t_s + 1;
    Dt = second(Timestamp(t_s,1)) -  second(Timestamp(t_s-1,1));
    s.x = [ meas(t_s,1); meas(t_s,2); (meas(t_s,1)-meas(t_s-1,1))*Dt; (meas(t_s,2)-meas(t_s-1,2))*Dt]; 
end


tic
s.x
kf = KalmanFilter_new(s);
t_1 = datevec(datestr(Timestamp(t_s,1)),'dd-mm-yyyy HH:MM:SS');
t_2 = datevec(datestr(Timestamp(end,1)),'dd-mm-yyyy HH:MM:SS');
v_t = t_1;
ind = t_s+1;
diff=etime(t_2,t_1);
for t=t_s+1:etime(t_2,t_1)
   v_t = datevec(addtodate(datenum(v_t),1,'second'))
   n_t = datevec(datestr(Timestamp(ind,1)),'dd-mm-yyyy HH:MM:SS')
   if (v_t==n_t)
       kf.s(end).A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
       kf.s(end).Q = [1/3, 0, 1/2, 0;  0, 1/3, 0, 1/2; 1/2, 0, 1, 0; 0, 1/2, 0, 1];
       kf.s(end+1) = kf.Predict(kf.s(end));
       kf.s(end).z = transpose(meas(ind,:)); % create a measurement
       ind = ind + 1
       kf.s(end) = kf.Update(kf.s(end));
   elseif v_t>n_t
       ind = ind + 1
   end
   z = transpose([kf.s(1:end-1).z]);
end
kf = kf.Smooth();
toc
figure
hold on
grid on
% plot measurement data:
hz=plot(z(:,1),z(:,2),'r.');
% plot a-posteriori state estimates:
x = transpose([kf.s(2:end).x]);
x_smooth = kf.x_smooth';
hk=plot(x(:,1),x(:,2),'b-');
x = xV';
hy=plot(x(:,1), x(:,2), 'g-');
%hy=quiver3(s.x(:,1),s.x(:,2),s.x(:,3),s.x(:,4));
%ht=plot(tru,'g-');
%legend([hz hk ht],'observations','Kalman output','true voltage',0)
title('Automobile Voltimeter Example')
