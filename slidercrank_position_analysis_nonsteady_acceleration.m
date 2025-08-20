% THIS CODE INCLUDES NON-CONSTANT VALUE OF CRANK ACCELERATION
% slidercrank_position_analysis_nonsteady_acceleration.m
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
N = 1000;
xB = zeros(2, N);
xC = zeros(2, N);
theta2 = zeros(1, N);
theta3 = zeros(1, N);
% allocate space for piston position
d = zeros(1, N);

% velocity analysis
omega2 = zeros(1, N); %not constant
omega3 = zeros(1, N);
ddot = zeros(1, N);

% acceleration analysis
alpha2 = zeros(1, N); %not constant
alpha3 = zeros(1, N);
d_doubledot = zeros(1, N);

% unsteady acceleration variables
A = 1; %pulse amplitude [rad/s/s]
T = 4; %period of pulse [sec]
t = zeros(1, N);
dt = 0.01; %increment

for i=1:N
    t(i) = (i-1)*dt;
    if (t(i) < T)
        %calculate crank angle, angular vel and accel for t < T
        theta2(i) = 0.5*A*t(i)^2 + T^2/(4*pi^2)*(cos(2*pi*t(i)/T) - 1);
        omega2(i) = A*t(i) - A*T/(2*pi)*sin((2*pi*t(i))/T);
        alpha2(i) = A*(1 - cos((2*pi*t(i))/T));
    else
        %calculate for t > T
        theta2(i) = A*T*t(i) - T^2*A/2;
        omega2(i) = A*T;
        alpha2(i) = 0;
    end

    theta3(i) = asin((c - a * sin(theta2(i))) / b);
    d(i) = a * cos(theta2(i)) + b * cos(theta3(i));

    % calculate unit vector
    [e1, n1] = UnitVector(0);
    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));

    % calculate position of point B and C
    xB(:, i) = FindPos(x0, a, e2);
    xC(:, i) = FindPos(xB(:, i), b, e3);

    % velocity analysis - solve for omega 3 and ddot
    A_mat = [b*n3, -e1];
    b_vec = -omega2(i)*a*n2;
    omega_vec = A_mat\b_vec;
    omega3(i) = omega_vec(1);
    ddot(i) = omega_vec(2);

    % acceleration analysis
    C_matrix = [b*n3, -e1];
    d_vector = -a*alpha2(i)*n2 + a*omega2(i)^2*e2 + b*omega3(i)^2*e3;
    alpha_vector = C_matrix\d_vector;
    alpha3(i) = alpha_vector(1);
    d_doubledot(i) = alpha_vector(2);
end

% validate result using numerical method
timestep = 2*pi/((N-1)*omega2(i));
d_doubledot_estimate = FiniteDiffMethod(ddot, timestep);

%convert XB and xC into mm
xB = xB * 1000;
xC = xC * 1000;
d = d * 1000;

%make three figures
fig1 = figure;
fig2 = figure;
fig3 = figure;
fig4 = figure;

%plot acceleration on fig4
figure(fig4);
%plot(t, alpha2);
plot(t, d_doubledot);
set(fig4, 'Position', [600, 100, 500, 300])
ylabel('Piston Acceleration')
title('Time vs Crank Acceleration');

%plot crank angle vs piston location on fig2
figure(fig2);
%plot(t, theta2);
plot(t, xC(1, :));
title('Time vs Crank Angle');
grid on
set(fig2, 'Position', [50, 500, 500, 300])

%plot velocity vs angle
figure(fig3);
%plot(t, omega2);
plot(t, ddot);
set(fig3, 'Position', [50, 100, 500, 300]);
grid on
title('Time vs Crank Velocity');

%now explicitly call fig1 so we can start plotting on it
figure(fig1);
set(fig1, 'Position', [600, 500, 500, 300])
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

while j <= N
    %avoid erratic behaviour if we click X to close plot
    if ~ishandle(fig1)
        break
    end

    set(h1, 'XData', [x0(1), xB(1,j), xC(1,j)], 'YData', [x0(2), xB(2, j), xC(2, j)]);
    set(h2, 'XData', [x0(1), xB(1,j), xC(1,j)], 'YData', [x0(2), xB(2, j), xC(2, j)]);
    drawnow

    j = j + 1;
    if loop == true & j >= N
        j = 1;
    end
end