% code save as from slidercrank_position_analysis.m
% for testing the force matrix developed for slider crank
% Afiq Rosli, 3/5/2025

% prepare workspace
clear variables; close all; clc;

% linkage dimensions
a = 0.040; % crank length (m)
b = 0.120; % connecting rod length (m)
c = 0.0;   % vertical slider offset (m)
p = 0.01;  % distance from point C to CM of piston

% inertial properties
m2 = 1.5; % crank mass (kg)
m3 = 0.08; % con rod mass (kg)
m4 = 0.16; % piston mass (kg)
I2 = 0.002; % crank moment of inertia (kg.m^2)
I3 = 0.0001; % con rod moment of inertia (kg.m^2) 
xbar2 = [0; 0]; % crank center of mass
xbar3 = [b/2; 0]; % con rod center of mass 
xbar4 = [p; 0]; % piston center of mass 

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

% velocity analysis
v0 = [0; 0]; %velocity of ground pin at A (origin)
omega2 = 10; 
omega3 = zeros(1, N);
ddot = zeros(1, N);

% acceleration analysis
a0 = [0; 0]; % acceleration of pin A
alpha2 = 0;
alpha3 = zeros(1, N);
d_doubledot = zeros(1, N);

% force analysis
FP = [-100; 0];
[vB, vC, v2, v3, v4] = deal(zeros(2, N));
[aB, aC, a2, a3, a4] = deal(zeros(2, N));
[FA, FB, FC, FD] = deal(zeros(2, N));
[M, T2] = deal(zeros(1, N));
[PExt, PKin] = deal(zeros(1, N));

U2 = eye(2);
Z2 = zeros(2);
Z21 = zeros(2, 1);
Z12 = zeros(1, 2);

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

    % velocity analysis - solve for omega 3 and ddot
    A_mat = [b*n3, -e1];
    b_vec = -omega2*a*n2;
    omega_vec = A_mat\b_vec;
    omega3(i) = omega_vec(1);
    ddot(i) = omega_vec(2);

    % acceleration analysis
    C_matrix = [b*n3, -e1];
    d_vector = -a*alpha2*n2 + a*omega2^2*e2 + b*omega3(i)^2*e3;
    alpha_vector = C_matrix\d_vector;
    alpha3(i) = alpha_vector(1);
    d_doubledot(i) = alpha_vector(2);

    aB(:, i) = FindAcc(a0, a, omega2, alpha2, e2, n2);
    aC(:, i) = FindAcc(aB(:, i), b, omega3(i), alpha3(i), e3, n3);

    % force analysis
    [eA2, nA2, LA2, s2A, s2B, ~] = LinkCG(a, 0, 0, xbar2, theta2(i));
    [eB3, nB3, LB3, s3B, s3C, ~] = LinkCG(b, 0, 0, xbar3, theta3(i));
    [eC4, nC4, LC4, s4C, s4D, ~] = LinkCG(p, 0, 0, xbar4, 0);

    a2(:, i) = FindAcc(a0, LA2, omega2, alpha2, eA2, nA2);
    a3(:, i) = FindAcc(aB(:, i), LB3, omega3(i), alpha3(i), eB3, nB3);
    a4(:, i) = FindAcc(aC(:, i), LC4, 0, 0, eC4, nC4);

    S_matrix = [U2    -U2    Z2    Z2    Z21    Z21;
                Z2     U2   -U2    Z2    Z21    Z21;
                Z2     Z2    U2   -U2    Z21    Z21;
              s2A'  -s2B'   Z12   Z12      0      1;
               Z12    s3B' -s3C'  Z12      0      0;
               Z12    Z12   s4C'  Z12      1      0;
               Z12    Z12   Z12   e1'      0      0];

    t_vector = [m2*a2(:, i);
                m3*a3(:, i);
                m4*a4(:, i) - FP;
                I2*alpha2;
                I3*alpha3(i);
                0;
                0];

    f_vector = S_matrix\t_vector;

    FA(:, i) = [f_vector(1); f_vector(2)];
    FB(:, i) = [f_vector(3); f_vector(4)];
    FC(:, i) = [f_vector(5); f_vector(6)];
    FD(:, i) = [f_vector(7); f_vector(8)];
    M(i) = f_vector(9);
    T2(i) = f_vector(10);

    % calculate velocity for energy method verification
    vB(:, i) = FindVel(v0, a, omega2, n2);
    vC(:, i) = FindVel(vB(:, i), b, omega3(i), n3);
    v2(:, i) = FindVel(v0, LA2, omega2, nA2);
    v3(:, i) = FindVel(vB(:, i), LB3, omega3(:, i), nB3);
    v4(:, i) = FindVel(vC(:, i), p, 0, eC4);

    %energy method verification
    PF = dot(FP, v4(:, i));
    PT = T2(i)*omega2;
    PExt(i) = PF + PT;

    P2 = InertialPower(m2, I2, v2(:, i), a2(:, i), omega2, alpha2);
    P3 = InertialPower(m3, I3, v3(:, i), a3(:, i), omega3(i), alpha3(i));
    P4 = InertialPower(m4, 0, v4(:, i), a4(:, i), 0, 0);
    PKin(i) = P2 + P3 + P4;
end

% validate result using numerical method
timestep = 2*pi/((N-1)*omega2);
d_doubledot_estimate = FiniteDiffMethod(ddot, timestep);

%convert XB and xC into mm
xB = xB * 1000;
xC = xC * 1000;
d = d * 1000;

layout = tiledlayout(3, 2);

%plot power
nexttile(4);
plot(theta2*180/pi, PKin);
hold on
plot(theta2*180/pi, PExt, '.');
xlim([0, 360])
ylabel('Power (W)');
xlabel(['Crank Angle(' char(176) ')']);
legend('Kinetic', 'External');

%plot torque
nexttile(6);
plot(theta2*180/pi, T2);
xlim([0 360]);
grid on
set(gca, 'xtick', 0:60:360);
set(gca, 'ytick', -5:1:5);
ylabel('Torque (Nm)');

%plot acceleration 
nexttile(5);
plot(theta2*180/pi, d_doubledot);
xlim([0 360]);
set(gca, 'xtick', 0:60:360);
ylabel('Piston Acceleration')
hold on
plot(theta2*180/pi, d_doubledot_estimate, '.')

%plot crank angle vs piston location 
nexttile(2);
plot(theta2 * 180 / pi, d);
title('Crank Angle vs Piston Location');
xlabel('Theta2 (degrees)');
ylabel('Piston location (mm)');
grid on
set(gca, 'xtick', 0:60:360);
xlim([0 360])

%plot velocity vs angle
nexttile(3);
plot(theta2(1, :) * 180/pi, ddot(1, :));
xlim([0 360]);
grid on
set(gca, 'xtick', 0:60:360);
xlabel('Theta2 (degrees)');
ylabel('Piston Velocity (mm/s)');
title('Crank Angle vs Piston Velocity');

%now explicitly call fig1 so we can start plotting on it
nexttile(1)
hold on;
plot(xB(1, :), xB(2, :));
plot(xC(1, :), xC(2, :));
axis equal
grid on
title('Paths of points B and C on the Slider-Crank')
xlabel('X-position (mm)');
ylabel('Y-position (mm)');
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
    if ~ishandle(layout)
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