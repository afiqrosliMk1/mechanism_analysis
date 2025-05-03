% slidercrank_position_analysis.m
% performs a position analysis on the slider-crank linkage and
% plots the piston position as a function of crank angle
% Afiq Rosli, 3/5/2025

% prepare workspace
clear variables; close all; clc;

% linkage dimensions
a = 0.040; % crank length (m)
b = 0.120; % connecting rod length (m)
c = 0.0;   % vertical slider offset (m)

% CHECK IF LINKAGE CAN BE ASSEMBLED BEFORE MOVING ON WITH THE CODE
if ~(a + abs(c) <= b)
    fprintf("CANNOT ASSEMBLE\n")
    return;
else
    fprintf("proceed\n")
end

% ground pins
x0 = [0; 0]; % ground pin at A (origin)

% number of times to perform calculation
N = 361;
xB = zeros(2, N);
xC = zeros(2, N);
theta2 = zeros(1, N);
theta3 = zeros(1, N);
% allocate space for piston position
d = zeros(1, N);

for i=1:N
    theta2(i) = (i - 1) * 2 * pi / (N - 1);
    theta3(i) = asin((c - a * sin(theta2(i))) / b);
    d(i) = a * cos(theta2(i)) + b * cos(theta3(i));

    % calculate unit vector
    [e1, n1] = UnitVector(0);
    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));

    % calculate position of point B and C
    xB(:, i) = FindPos(x0, a, e2);
    xC(:, i) = FindPos(xB(:, i), b, e3);
end

%convert XB and xC into mm
xB = xB * 1000;
xC = xC * 1000;
d = d * 1000;

%make two figures
fig1 = figure;
fig2 = figure;
%plot crank angle vs piston location on fig2
plot(theta2 * 180 / pi, d);
title('Crank Angle vs Piston Location');
xlabel('Theta2 (degrees)');
ylabel('Piston location (mm)');
grid on
set(gca, 'xtick', 0:60:360);
set(fig2, 'Position', [50, 400, 600, 400])
xlim([0 360])
%now explicitly call fig1 so we can start plotting on it
figure(fig1);
set(fig1, 'Position', [700, 400, 600, 400])
hold on;
plot(xB(1, :), xB(2, :));
plot(xC(1, :), xC(2, :));
axis equal
grid on
title('Paths of points B and C on the Slider-Crank')
xlabel('X-position (mm)');
ylabel('Y-position (mm)');
%xlim([-0.05, 0.2]);
%ylim([-0.1, 0.1]);
legend('Point B', 'Point C', Location='eastoutside')

%flag to keep looping
loop = true;
j = 1;

%initialise linkage
h1 = plot([x0(1), xB(1,1), xC(1,1)], [x0(2), xB(2, 1), xC(2, 1)], LineWidth=1, Color='b', HandleVisibility='off');

%initialise joint marker
h2 = plot([x0(1), xB(1,1), xC(1,1)], [x0(2), xB(2, 1), xC(2, 1)], 'o', 'MarkerFaceColor','k', HandleVisibility='off');

while j <= N;
    %avoid erratic behaviour if we click X to close plot
    if ~ishandle(fig1)
        break
    end

    set(h1, 'XData', [x0(1), xB(1,j), xC(1,j)], 'YData', [x0(2), xB(2, j), xC(2, j)]);
    set(h2, 'XData', [x0(1), xB(1,j), xC(1,j)], 'YData', [x0(2), xB(2, j), xC(2, j)]);
    drawnow

    j = j + 1;
    if loop == true & j >= N;
        j = 1;
    end
end