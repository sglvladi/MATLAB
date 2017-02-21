lat=36.0762253*10;
lon = 24.87468932*10;
alt = 0;

a = 6378137;
b = 8.1819190842622e-2;

% intermediate calculation
% (prime vertical radius of curvature)
N = a / sqrt(1 - b^2 * sin(lat)^2);

% results:
x = (N+alt) * cos(lat) * cos(lon)
y = (N+alt) * cos(lat) * sin(lon)
z = ((1-b^2) * N + alt) * sin(lat)
