% threebar_position_analysis.m
% Conducts a position analysis on the threebar crank-slider linkage
% by Afiq Rosli, 17/11/2024

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 0.150; % crank length (m)
d = 0.400; % length between ground pins (m)
p = 0.300; % slider length (m)

% Ground pins. convention: variable starting with "x" is to store
% coordinate of points
x0 = [0; 0]; % point A (the origin) - vector
xD = [d; 0]; % point D - vector

% number of times to perform position calculations. 0-360 increment 1, will
% give you 361 calculations
N = 100;

theta2 = zeros(1, N); % allocate space for crank angle
theta3 = zeros(1, N); % allocate space for slider angle
b = zeros(1, N); % allocate space for slider length

xB = zeros(2, N);
xP = zeros(2, N);
xQ = zeros(2, N); %point where door is pushed

T = 2; %time it takes to open door
lambda = 2*pi/T;
delta_time = T/(N-1);
t = zeros(1, N);

% velocity analysis
omega2 = zeros(1, N) ; % rad/s
omega3 = zeros(1, N);
bdot = zeros(1, N);
v0 = [0; 0]; % velocity of point A

vB = zeros(2, N);
vP = zeros(2, N);

% acceleration analysis variables
alpha2 = zeros(1, N);
alpha3 = zeros(1, N);
b_doubledot = zeros(1, N);

a0 = [0; 0];
aB = zeros(2, N);
aP = zeros(2, N);

% force analysis variables
W = 1.0; %door width
L = 0.05; %door thickness
H = 2.0; %door height
m2 = 60; 
m3 = 0.0; %slider mass neglected. too small.
I2 = m2/12*(W^2 + L^2);
I3 = 0.0;
xbar2 = [W/2; 0];
xbar3 = [p/2; 0];

k = 1000; %spring constant [N/m]
b0 = d - a; %unstretched spring length

U2 = eye(2);
Z2 = zeros(2, 2);
Z12 = zeros(1, 2);

[FA, FB, FD, FQ] = deal(zeros(2, N));

[v2, v3, vQ] = deal(zeros(2, N));
[a2, a3] = deal(zeros(2, N));
[PExt, P_Kin] = deal(zeros(1, N));

% Main Loop
for i = 1:N

    t(i) = (i-1) * delta_time; %minus 1 force the first t to be zero
    if t(i) <= T
        theta2(i) = 0.25*(lambda*t(i) - sin(lambda*t(i)));
        omega2(i) = 0.25*lambda*(1 - cos(lambda*t(i)));
        alpha2(i) = 0.25*lambda^2*sin(lambda*t(i));
    elseif t(i) > T
        fprintf("t exceeding specified period T. function to handle this not included yet\n");
    end
    

    theta3(i) = atan2(-a * sin(theta2(i)) , d - a * cos(theta2(i)));
    b(i) = (d - a * cos(theta2(i))) / cos(theta3(i));

    %calculate unit vectors
    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));

    %solve for positions of points B and P on the linkage
    xB(:, i) = FindPos(x0, a, e2);
    xP(:, i) = FindPos(xB(:, i), p, e3); 
    xQ(:, i) = FindPos(x0, W, e2); %the point where door is pushed

    %velocity analysis - solve for omega3 and bdot
    A_mat = [b(i)*n3, e3];
    b_vec = -omega2(i)*a*n2;
    omega_vec = A_mat\b_vec;
    omega3(i) = omega_vec(1);
    bdot(i) = omega_vec(2);

    %velocity analysis - solve for point velocity
    vB(:, i) = FindVel(v0, a, omega2(i), n2);
    vP(:, i) = FindVel(vB(:, i), p, omega3(i), n3);

    %acceleration analysis - solve for alpha3 and b_doubledot
    C_matrix = [b(i)*n3, e3];
    d_vector = -a*alpha2(i)*n2 + a*omega2(i)^2*e2 - 2*bdot(i)*omega3(i)*n3 + b(i)*omega3(i)^2*e3;
    alpha_vector = C_matrix\d_vector;
    alpha3(i) = alpha_vector(1);
    b_doubledot(i) = alpha_vector(2);

    %acceleration analysis - find point acceleration
    aB(:, i) = FindAcc(a0, a, omega2(i), alpha2(i), e2, n2);
    aP(:, i) = FindAcc(aB(:, i), p, omega3(i), alpha3(i), e3, n3);

    %force analysis
    [eA2, nA2, LA2, s2A, s2B, s2Q] = LinkCG(a, W, 0, xbar2, theta2(i) );
    [eB3, nB3, LB3, s3B, s3D, s3P] = LinkCG(b(i), p, 0, xbar3, theta3(i));

    v2(:, i) = FindVel(v0, LA2, omega2(i), nA2);
    v3(:, i) = FindVel(vB(:, i), LB3, omega3(i), nB3);
    vQ(:, i) = FindVel(v0, W, omega2(i), n2);

    a2(:,i) = FindAcc(a0, LA2, omega2(i), alpha2(i), eA2, nA2);
    a3(:, i) = FindAcc(aB(:, i), LB3, omega3(i), alpha3(i), eB3, nB3);

    Fk = k*(b(i) - b0)*e3;

    S_matrix = [U2    -U2    Z2    U2;
                Z2     U2   -U2    Z2;
                s2A'  -s2B'  Z12   s2Q';
                Z12    s3B' -s3D'  Z12;
                Z12    Z12   e3'   Z12;
                Z12    Z12   Z12   e2'];

    t_vector = [m2*a2(:, i);
                m3*a3(:, i) - Fk;
                I2*alpha2(i);
                I3*alpha3(i);
                0;
                0];

    f_vector = S_matrix\t_vector;

    FA(:, i) = [f_vector(1); f_vector(2)];
    FB(:, i) = [f_vector(3); f_vector(4)];
    FD(:, i) = [f_vector(5); f_vector(6)];
    FQ(:, i) = [f_vector(7); f_vector(8)];

    %verify result using energy method
    PQ = dot(FQ(:, i), vQ(:, i));
    PK = dot(Fk, vB(:, i));
    PExt(i) = PQ + PK;

    P2 = InertialPower(m2, I2, v2(:, i), a2(:, i), omega2(i), alpha2(i));
    P3 = InertialPower(m3, I3, v3(:, i), a3(:, i), omega3(i), alpha3(i));
    P_Kin(i) = P2 + P3;

end

%loop to estimate velocity and acceleration using numerical method
%---purely just to validate result of analytical vs numerical method-----
%velocity
%dt = 2*pi/((N-1)*omega2); 
dt = delta_time; %we are using delta_time value for dt
omega3_estimate = FiniteDiffMethod(theta3, dt);
bdot_estimate = FiniteDiffMethod(b, dt);
vPx_estimate = FiniteDiffMethod(xP(1, :), dt);

%acceleration
aPx_estimate = FiniteDiffMethod(vP(1, :), dt);
%-----------

loop = true;

fig = figure('Units','normalized','Position',[0.2 0.1 0.35 0.8]);
layout = tiledlayout(2, 2);
nexttile
grid on
%plot path B and P and Q
plot(xB(1, :), xB(2, :));
hold on
plot(xP(1, :), xP(2, :));
hold on
plot(xQ(1, :), xQ(2, :), '--');
axis equal
xlim([0, 1.1]);
ylim([0, 1.1]);

title('Paths of points B and P on the Threebar Linkage')
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point B', 'Point P', 'Location', 'eastoutside')

%initialise linkage plot 
x = [x0(1), xB(1, 1), xP(1, 1)];
y = [x0(2), xB(2, 1), xP(2, 1)];
%handle visibility=off prevents h1 from being added to the legend
h1 = plot(x, y, LineWidth=2, Color='k', HandleVisibility='off');

%plot the door
h3 = plot([xB(1, 1), xQ(1, 1)], [xB(1, 1), xQ(2, 1)], HandleVisibility='off');

%initialise pin joint
h2 = plot([x0(1), xD(1), xB(1, 1), xP(1, 1), xQ(1, 1)], [x0(2), xD(2), ...
        xB(2, 1), xP(2, 1), xQ(2, 1)], 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'b', HandleVisibility='off');

%initialise text object
tA = text(x0(1) - 0.015, x0(2), 'A', HorizontalAlignment='center');
tB = text(xB(1, 1) - 0.015, xB(2, 1) - 0.015, 'B', HorizontalAlignment='center', VerticalAlignment='middle');
tD = text(xD(1) - 0.015, xD(2), 'D', HorizontalAlignment='center');
tP = text(xP(1, 1) - 0.015, xP(2, 1), 'P', HorizontalAlignment='center');
tQ = text(xQ(1, 1) - 0.015, xQ(2, 1), 'Q', HorizontalAlignment='center');

%plot for velocity
nexttile
plot(theta2*180/pi, vP(1, :));
xlabel('angle (degree)')
xlim("auto");
set(gca, 'xtick', 0:5:90);

%----lines below is for plotting estimated value from numerical method
hold on
plot(theta2*180/pi, vPx_estimate(1, :), '.');
legend('vPx', 'vPx\_estimate', 'Location', 'best');

%--plot for acceleration
nexttile
plot(theta2*180/pi, aP(1, :));
xlim("auto");
set(gca, 'xtick', 0:5:90);
hold on
plot(theta2*180/pi, aPx_estimate(1, :), '.');
legend('aPx', 'aPx\_estimate', 'Location', 'best');

%plot for force
nexttile
plot(t, FQ(1, :));
hold on
plot(t, FQ(2, :));
hold on
plot(t, sqrt(sum(FQ.^2, 1)));
xlabel('time (s)');
ylabel('Force (N)');
legend('FQ_x', 'FQ_y', 'FQ');

% plot for verify force result using energy method
fig2 = figure; 
plot(t, P_Kin);
hold on
plot(t, PExt, '.');
xlabel('time (s)');
ylabel('Power (W)');
title('Kinetic vs External Power');
legend('Kinetic', 'External');

j = 1;
while j <= N
    %avoid erratic behaviour if we click X to close plot
    if ~ishandle(fig)
        break
    end

    %update linkage position
    set(h1, "XData", [x0(1), xB(1, j), xP(1, j)], "YData", [x0(2), xB(2, j), xP(2, j)]);
    set(h2, "XData", [x0(1), xD(1), xB(1, j), xP(1, j), xQ(1, j)], "YData", [x0(2), xD(2), xB(2, j), xP(2, j), xQ(2, j)]);
    set(h3, "XData", [xB(1, j), xQ(1, j)], "YData", [xB(2, j), xQ(2, j)]);

    %update labels position
    set(tB, 'Position', [xB(1, j) - 0.015, xB(2, j) - 0.015, 0])
    set(tP, 'Position', [xP(1, j) - 0.015, xP(2, j) - 0.015, 0])
    set(tQ, 'Position', [xQ(1, j) - 0.015, xQ(2, j) - 0.015, 0]);

    drawnow

    j = j + 1;
    if j == N & loop == true
        j = 1;
    end

end