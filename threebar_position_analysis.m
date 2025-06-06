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

loop = true;

fig = figure;
hold on
grid on
%plot path B and P
plot(xB(1, :), xB(2, :));
plot(xP(1, :), xP(2, :));
axis equal
xlim([-0.4, 0.4]);
ylim([-0.15, 0.15]);

title('Paths of points B and P on the Threebar Linkage')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point B', 'Point P', 'Location', 'eastoutside')

%initialise linkage plot 
x = [x0(1), xB(1, 1), xP(1, 1)];
y = [x0(2), xB(2, 1), xP(2, 1)];
%handle visibility=off prevents h1 from being added to the legend
h1 = plot(x, y, LineWidth=2, Color='k', HandleVisibility='off');

%initialise pin joint
h2 = plot([x0(1), xD(1), xB(1, 1), xP(1, 1)], [x0(2), xD(2), xB(2, 1),...
        xP(2, 1)], 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b', HandleVisibility='off');

%initialise text object
tA = text(x0(1) - 0.015, x0(2), 'A', HorizontalAlignment='center');
tB = text(xB(1, 1) - 0.015, xB(2, 1) - 0.015, 'B', HorizontalAlignment='center', VerticalAlignment='middle');
tD = text(xD(1) - 0.015, xD(2), 'D', HorizontalAlignment='center');
tP = text(xP(1, 1) - 0.015, xP(2, 1), 'P', HorizontalAlignment='center');

j = 1;
while j <= N
    %avoid erratic behaviour if we click X to close plot
    if ~ishandle(fig)
        break
    end

    %update linkage position
    set(h1, "XData", [x0(1), xB(1, j), xP(1, j)], "YData", [x0(2), xB(2, j), xP(2, j)]);
    set(h2, "XData", [x0(1), xD(1), xB(1, j), xP(1, j)], "YData", [x0(2), xD(2), xB(2, j), xP(2, j)]);

    %update labels position
    set(tB, 'Position', [xB(1, j) - 0.015, xB(2, j) - 0.015, 0])
    set(tP, 'Position', [xP(1, j) - 0.015, xP(2, j) - 0.015, 0])

    drawnow

    j = j + 1;
    if j == N & loop == true
        j = 1;
    end


end