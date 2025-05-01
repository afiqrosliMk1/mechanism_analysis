% threebar_position_analysis.m
% Conducts a position analysis on the threebar crank-slider linkage
% by Afiq Rosli, 17/11/2024

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 0.100; % crank length (m)
d = 0.150; % length between ground pins (m)
p = 0.300; % slider length (m)

% Ground pins. convention: variable starting with "x" is to store
% coordinate of points
x0 = [0; 0]; % point A (the origin) - vector
xD = [d; 0]; % point D - vector

% number of times to perform position calculations. 0-360 increment 1, will
% give you 361 calculations
N = 361; 

theta2 = zeros(1, N); % allocate space for crank angle
theta3 = zeros(1, N); % allocate space for crank angle
b = zeros(1, N); % allocate space for slider length

xB = zeros(2, N);
xP = zeros(2, N);

% Main Loop
for i = 1:N
    %(2*pi)/(N-1) gives you angle per step. 
    theta2(i) = (i - 1) * (2 * pi) / (N - 1);
    theta3(i) = atan2(-a * sin(theta2(i)) , d - a * cos(theta2(i)));
    %theta3(i) = atan2(a * sin(theta2(i)) , a * cos(theta2(i)) - d);
    b(i) = (d - a * cos(theta2(i))) / cos(theta3(i));

    %calculate unit vectors
    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));

    %solve for positions of points B and P on the linkage
    xB(:, i) = FindPos(x0, a, e2);
    xP(:, i) = FindPos(xB(:, i), p, e3); 
end

plot(xB(1, :), xB(2, :))
hold on
plot(xP(1, :), xP(2, :))
axis equal
grid on

title('Paths of points B and P on the Threebar Linkage')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point B', 'Point P', 'Location', 'southeast')

%draw linkage
x = [x0(1), xB(1, 80), xP(1, 80)];
y = [x0(2), xB(2, 80), xP(2, 80)];
plot(x, y, LineWidth=2, Color='k');

%draw linkage joint
plot([x0(1), xD(1), xB(1, 80), xP(1, 80)], [x0(2), xD(2), xB(2, 80),...
    xP(2, 80)], 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'b');

%draw labels on each pin joint
text(x0(1) - 0.015, x0(2), 'A', HorizontalAlignment='center')
text(xB(1, 80) - 0.015, xB(2, 80), 'B', HorizontalAlignment='center')
text(xD(1) - 0.015, xD(2), 'D', HorizontalAlignment='center')
text(xP(1, 80) - 0.015, xP(2, 80), 'P', HorizontalAlignment='center')
