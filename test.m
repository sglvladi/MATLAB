[m,n] = size(lat);
measurements = zeros(m,3);
%for i=1:m 
[measurements(:,1), measurements(:,2),measurements(:,3)] = lla2ecef(lat(:,1)*10,long(:,1)*10,0);
%end
    