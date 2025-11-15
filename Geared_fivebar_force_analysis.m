% Gear_fivebar_force_analysis.m
% saved from Gear_fivebar.m
% conducts force analysis
% Afiq Rosli, 09/11/2025

% prepare workspace
clear variables; close all; clc;
a = 0.120;               % crank length (m)
b = 0.250;               % coupler 1 length (m)
c = 0.250;               % coupler 2 length (m)
d = 0.180;               % ground link length (m)
u = 0.120;               % length of link on gear 2 (m)
N1 = 48;                 % number of teeth on gear 1
N2 = 24;                 % number of teeth on gear 2
rho = N1/N2;             % gear ratio
r1 = d*rho/(rho + 1);    % radius of gear 1
r2 = d/(rho + 1);        % radius of gear 2
phi = 0;                 % offset angle between gears
gamma3 = 20 * pi / 180;  % angle to point P on coupler 1
gamma4 = -20 * pi / 180; % angle to point Q on coupler 2
p = 0.200;               % distance to point P on coupler 1
q = 0.200;               % distance to point Q on coupler 2
xbar2 = [a/2; 0];
xbar3 = [b/2; 0];
xbar4 = [c/2; 0];
xbar5 = [u/2; 0];
pressure_angle = 20*pi/180; % pressure angle 20degree
m2 = 0.100; % mass of gear 1 (kg)
m3 = 0.200; % mass of coupler 1 (kg)
m4 = 0.200; % mass of coupler 2 (kg)
m5 = 0.100; % mass of gear 2 (kg)
I2 = 0.0001; % moment of inertia of gear 1 about cg
I3 = 0.0002; % moment of inertia of coupler 1 about cg
I4 = 0.0002; % moment of inertia of coupler 2 about cg
I5 = 0.0001; % moment of inertia of gear 2 about cg

I2A = I2 + m2*dot(xbar2, xbar2); % moment of inertia of gear 1 about pin A (kgm^2)
I5D = I5 + m5*dot(xbar5, xbar5); % moment of inertia of gear 2 about pin D

% % check for grashoff condition. terminate if non-grashoff
% [S, L, P, Q] = slpq(a , b, c, d);
% if S + L < P + Q 
%     fprintf("grashoff\n")
% else
%     fprintf("non-grashoff\n")
%     return
% end

% initialise variables
loop = true;
N = 361;
theta2 = zeros(1, N);
theta3 = zeros(1, N);
theta4 = zeros(1, N);
delta = 0;
r = 0;
s = 0;
f = 0;
g = 0;
h = 0;
x0 = [0; 0];
xD = [d; 0];
[xB, xC, xP, xE, xQ] = deal(zeros(2, N));

% variables for velocity analysis
omega2 = 10; %rad/s
omega3 = zeros(1, N);
omega4 = zeros(1, N);
omega5 = -rho*omega2;

v0 = [0; 0]; %velocity at crank pin
[vB, vC, vP, vQ] = deal(zeros(2, N));
[v2, v3, v4, v5] = deal(zeros(2, N));

% variables for acceleration analysis
alpha2 = 0; %rad/s/s
alpha5 = -rho*alpha2;
alpha3 = zeros(1, N);
alpha4 = zeros(1, N);

a0 = [0; 0]; %acceleration of origin
[aB, aC, aP] = deal(zeros(2, N));
[a2, a3, a4, a5] = deal(zeros(2, N));

% variables for force analysis
U2 = eye(2);
Z2 = zeros(2);
Z12 = zeros(1, 2);
Z21 = zeros(2, 1);

FQ = [0; -100];

[FA, FB, FC, FD, FE, Fg] = deal(zeros(2, N));
T2 = zeros(1, N);

[P_Ext, P_Kin] = deal(zeros(1, N));


% main loop
for i = 1:N 
     theta2(i) = (i - 1) * (2*pi) / (N - 1); 
     theta5 = -rho * theta2(i) + phi;  % angle of second gear
     [e2, n2] = UnitVector(theta2(i));
     [e5, n5] = UnitVector(theta5);
     xC(:, i) = FindPos(xD, u, e5);

     dprime = sqrt(xC(1, i)^2  + xC(2, i)^2);
     beta = atan2(xC(2, i), xC(1, i));

     r = dprime - a * cos(theta2(i) - beta);
     s = a * sin(theta2(i) - beta);

     f = sqrt(r^2 + s^2);

     % open configuration
     delta = acos((b^2 + c^2 - f^2)/(2*b*c));

     g = b - c * cos(delta);
     h = c * sin(delta);

     theta3(i) = atan2((h*r - g*s), (g*r + h*s)) + beta;

     theta4(i) = theta3(i) + delta;

     
     [e3, n3] = UnitVector(theta3(i));
     [e4, n4] = UnitVector(theta4(i));
     [eBP, nBP] = UnitVector(theta3(i) + gamma3);
     [eCQ, nCQ] = UnitVector(theta4(i) + gamma4);

     xB(:, i) = FindPos(x0, a, e2);
     xE(:, i) = FindPos(xB(:, i), b , e3);
     xP(:, i) = FindPos(xB(:, i), p, eBP);
     xQ(:, i) = FindPos(xC(:, i), q, eCQ);

     % velocity analysis
     A_matrix = [b*n3, -c*n4];
     b_vector = -a*omega2*n2 + u*omega5*n5;
     omega_matrix = A_matrix \ b_vector;
     omega3(i) = omega_matrix(1);
     omega4(i) = omega_matrix(2);

     vB(:, i) = FindVel(v0, a, omega2, n2);
     vP(:, i) = FindVel(vB(:, i), p, omega3(i), nBP);
     vC(:, i) = FindVel(v0, u, omega5, n5);
     vQ(:, i) = FindVel(vC(:, i), q, omega4(i), nCQ);

     % acceleration analysis
     C_matrix = [b*n3, -c*n4];
     d_vector = -a*alpha2*n2 + a*omega2^2*e2 + b*omega3(i)^2*e3 - ...
         c*omega4(i)^2*e4 +u*alpha5*n5 - u*omega5^2*e5;
     alpha_matrix = C_matrix\d_vector;
     alpha3(i) = alpha_matrix(1);
     alpha4(i) = alpha_matrix(2);

     aB(:, i) = FindAcc(a0, a, omega2, alpha2, e2, n2);
     aP(:, i) = FindAcc(aB(:, i), p, omega3(i), alpha3(i), eBP, nBP);
     aC(:, i) = FindAcc(a0, u, omega5, alpha5, e5, n5);

     % force analysis
     [eA2, nA2, LA2, s2A, s2B, ~] = LinkCG(a, 0, 0, xbar2, theta2(i));
     [eB3, nB3, LB3, s3B, s3E, s3P] = LinkCG(b, p, gamma3, xbar3, theta3(i));
     [eC4, nC4, LC4, s4C, s4E, s4Q] = LinkCG(c, q, gamma4, xbar4, theta4(i));
     [eD5, nD5, LD5, s5D, s5C, ~] = LinkCG(u, 0, 0, xbar5, theta5);

     v2(:, i) = FindVel(v0, LA2, omega2, nA2);
     v3(:, i) = FindVel(vB(:, i), LB3, omega3(i), nB3);
     v4(:, i) = FindVel(vC(:, i), LC4, omega4(i), nC4);
     v5(:, i) = FindVel(v0, LD5, omega5, nD5);


     sAB = a*n2;
     sDC = u*n5;

     a2(:, i) = FindAcc(a0, LA2, omega2, alpha2, eA2, nA2);
     a3(:, i) = FindAcc(aB(:, i), LB3, omega3(i), alpha3(i), eB3, nB3);
     a4(:, i) = FindAcc(aC(:, i), LC4, omega4(i), alpha4(i), eC4, nC4);
     a5(:, i) = FindAcc(a0, LD5, omega5, alpha5, eD5, nD5);

     S_matrix = [U2     Z2   -U2    Z21    Z21;
                 Z2    -U2    U2    Z21    Z21;
                 -sAB'  Z12   Z12    r1      1;
                 s3B'  Z12  -s3E'    0      0;
                 Z12   -s4C' s4E'    0      0;
                 Z12   sDC'  Z12    r2      0;];

     t_vector = [m3*a3(:, i);
                 m4*a4(:, i) - FQ;
                 I2A*alpha2;
                 I3*alpha3(i);
                 I4*alpha4(i) - dot(s4Q, FQ);
                 I5D*alpha5];

     f_vector = S_matrix\t_vector;

     FB(:, i) = [f_vector(1), f_vector(2)];
     FC(:, i) = [f_vector(3), f_vector(4)];
     FE(:, i) = [f_vector(5), f_vector(6)];
     Fg(2, i) = f_vector(7);
     T2(i) = f_vector(8);

     Fg(1, i) = -abs(Fg(2, i)*tan(pressure_angle));

     FA(:, i) = FB(:, i) + m2*a2(:, i) - Fg(:, i);
     FD(:, i) = FC(:, i) - m5*a5(:, i) - Fg(:, i);

     % Result validation using energy method    
     P_FQ = dot(FQ, vQ(:, i));
     P_T2 = T2(i)*omega2;
     P_Ext(i) = P_FQ + P_T2;

     P2 = InertialPower(m2, I2, v2(:, i), a2(:, i), omega2, alpha2);
     P3 = InertialPower(m3, I3, v3(:, i), a3(:, i), omega3(i), alpha3(i));
     P4 = InertialPower(m4, I4, v4(:, i), a4(:, i), omega4(i), alpha4(i));
     P5 = InertialPower(m5, I5, v5(:, i), a5(:, i), omega5, alpha5);
     P_Kin(i) = P2 + P3 + P4 + P5;
end

%numerical estimation of velocity and acceleration
timestep = 2*pi/((N-1)*omega2);
vPx_estimate = FiniteDiffMethod(xP(1, :), timestep);
vPy_estimate = FiniteDiffMethod(xP(2, :), timestep);

aPx_estimate = FiniteDiffMethod(vP(1, :), timestep);
aPy_estimate = FiniteDiffMethod(vP(2, :), timestep);

%create a tiled plot
layout = tiledlayout(3, 2);

%plot velocity or acceleration
nexttile(2)
plot(theta2*180/pi, aP(1, :));
hold on
plot(theta2*180/pi, aP(2, :));
set(gca, 'xtick', 0:60:360);
xlim([0, 360]);
%plot the estimated value for validating result
plot(theta2*180/pi, aPx_estimate, '.');
plot(theta2*180/pi, aPy_estimate, '.');
legend('aPx', 'aPy', 'aPx\_estimate', 'aPy\_estimate', 'Location', 'best');
ylabel('Acceleration (m/s^2)');
xlabel('Crank Angle (degree)');

%plot gear force
nexttile(4);
plot(theta2*180/pi, Fg(1, :));
hold on
plot(theta2*180/pi, Fg(2, :));
xlim([0 360]);
set(gca, 'xtick', 0:60:360);
grid on
ylabel('Gear Force (N)');
xlabel('Crank Angle (degree)');
legend('Fgx', 'Fgy', Location='best');

%plot torque
nexttile(5);
plot(theta2*180/pi, T2);
xlim([0 360]);
set(gca, 'xtick', 0:60:360);
grid on
ylabel('Torque (Nm)');
xlabel('Crank Angle (degree)');
legend('T2');

%plot power
nexttile(3)
plot(theta2*180/pi, P_Ext);
hold on
plot(theta2*180/pi, P_Kin, '.');
xlim([0 360]);
set(gca, 'xtick', 0:60:360);
grid on
ylabel('Power (W)');
xlabel('Crank Angle (degree)');
legend('P\_Ext', 'P\_Kin', Location='best');

% plot path
nexttile(1)
% plot(xB(1, :), xB(2, :));
plot(xQ(1, :), xQ(2, :));
hold on
plot(xP(1, :), xP(2, :));
%plot(xE(1, :), xE(2, :));

axis equal
xlim([-0.30, 0.40]);
ylim([-0.15, 0.4]);

%initialise linkage plot
h1 = plot([x0(1), xB(1, 1), xE(1, 1)], [x0(2), xB(2, 1), xE(2, 1)]);
h2 = plot([xB(1, 1), xP(1, 1), xE(1, 1)], [xB(2, 1), xP(2, 1), xE(2, 1)]);
h3 = plot([xD(1), xC(1, 1), xE(1, 1)], [xD(2), xC(2, 1), xE(2, 1)]);
h4 = patch([xB(1, 1), xP(1, 1), xE(1, 1)], [xB(2, 1), xP(2, 1), xE(2, 1)], [0.5, 0.5, 1]);
h6 = plot([xC(1, 1), xQ(1, 1), xE(1, 1)], [xC(2, 1), xQ(2, 1), xE(2, 1)]);
h7 = patch([xC(1, 1), xQ(1, 1), xE(1, 1)], [xC(2, 1), xQ(2, 1), xE(2, 1)], [1, 0.5, 0.5]);

%initialise joint marker
h5 = plot([x0(1), xB(1, 1), xE(1, 1), xD(1), xC(1, 1)], ...
     [x0(2), xB(2, 1), xE(2, 1), xD(2), xC(2, 1)], 'o', 'MarkerFaceColor', 'k');

%initialise label
t1 = text(x0(1) + 0.005, x0(2), 0, "A");
t2 = text(xB(1, 1) + 0.005, xB(2, 1), 0, "B");
t3 = text(xC(1, 1) + 0.005, xC(2, 1), 0, "C");
t4 = text(xE(1, 1) + 0.005, xE(2, 1), 0, "E");
t5 = text(xD(1) + 0.005, xD(2), 0, "D");
t6 = text(xP(1) + 0.005, xP(2), 0, "P");
t7 = text(xQ(1) + 0.005, xQ(2), 0, "Q");


% animation loop
j = 1;
while j <= N
    if ~ishandle(layout)
        break
    end
    set(h1, 'XData', [x0(1), xB(1, j), xE(1, j)], 'YData', [x0(2), xB(2, j), xE(2, j)]);
    set(h2, 'XData', [xB(1, j), xP(1, j), xE(1, j)], 'YData', [xB(2, j), xP(2, j), xE(2, j)]);
    set(h3, 'XData', [xD(1), xC(1, j), xE(1, j)], 'YData', [xD(2), xC(2, j), xE(2, j)]);
    set(h6, 'XData', [xC(1, j), xQ(1, j), xE(1, j)], 'YData', [xC(2, j), xQ(2, j), xE(2, j)]);

    %update joint marker
    set(h5, 'XData', [x0(1), xB(1, j), xE(1, j), xD(1), xC(1, j)], ...
    'YData', [x0(2), xB(2, j), xE(2, j), xD(2), xC(2, j)])

    %update coupler patch
    set(h4, 'XData', [xB(1, j), xP(1, j), xE(1, j)], 'YData', [xB(2, j), xP(2, j), xE(2, j)]);
    set(h7, 'XData', [xC(1, j), xQ(1, j), xE(1, j)], 'YData', [xC(2, j), xQ(2, j), xE(2, j)]);
    drawnow

    %check all text object still exist before update label position
    if all(isgraphics([t2, t3, t4, t6, t7]))
        set(t2, 'Position', [xB(1, j) + 0.005, xB(2, j)]);
        set(t3, 'Position', [xC(1, j) + 0.005, xC(2, j)]);
        set(t4, 'Position', [xE(1, j) + 0.005, xE(2, j)]);
        set(t6, 'Position', [xP(1, j) + 0.005, xP(2, j)]);
        set(t7, 'Position', [xQ(1, j) + 0.005, xQ(2, j)]);
    end

    j = j + 1;
    if j == N & loop == true
        j = 1;
    end
end

%function slpq 
function [S, L, P, Q] = slpq(a, b, c, d)

    % initialise vector A to store linkage lengths passed to the function
    A = zeros(1, 4);
    A(1) = a;
    A(2) = b;
    A(3) = c;
    A(4) = d;
    
    % get min value and assign to S along with the index
    [S, idx] = min(A);
    % delete element at the given index
    A(idx) = [];
    [L, idx] = max(A);
    A(idx) =[];
    % assign P and Q to the rest of the values left
    P = A(1);
    Q = A(2);

end