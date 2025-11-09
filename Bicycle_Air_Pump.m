% Bicycle_Air_Pump.m
% conducts force analysis on air pump inverted slider crank
% to get the air pump configureation, delta is 180deg, rocker length=0.
% Afiq Rosli,2/11/2025

% prepare workspace
clear variables; close all; clc;
a = 0.070;               % crank length (m)
c = 0.0;               % rocker length (m)
d = 0.250;              % ground link length (m)
p = 0.145;               % length from B to P (m)
q = 0.200;               % length from A to Q (m)
delta = 180 * pi/ 180;  % angle between slider and rocker (radian)
gamma3 = 0;             % angle to point P on piston
H = 0.025;         %distance from point D to top of cylinder
h = c * sin(delta); % h is constant. only calculate it once
eps = 0.000001; % tiny number added to theta2 to keep its bounds

% if grashoff condition met, grashoff = true
if ((a + d) < c * sin(delta))
    fprintf("linkage cannot be assembled\n")
else
    if abs(a - d) >= c * sin(delta)
        grashoff = true;
        fprintf("linkage is grashoff\n")
        % theta2_min = 0;
        % theta2_max = 2*pi;
    else
        grashoff = false;
        fprintf("non grashoff\n")
        % theta2_min = acos((a^2 + d^2 - (c*sin(delta))^2) / (2*a*d))  + eps;
        % theta2_max = 2*pi - theta2_min;
    end
end

% initialise variables
loop = true;
N = 361;
b = zeros(1, N);
theta2 = zeros(1, N);
theta3 = zeros(1, N);
theta4 = zeros(1, N);
r = 0;
s = 0;
f = 0;
g = 0;
x0 = [0; 0]; % ground pin at A (origin)
xD = [d; 0]; % ground pin at D
xB = zeros(2, N); 
xC = zeros(2, N);
xP = zeros(2, N);
xQ = zeros(2, N);

%variable for velocity analysis
omega2 = 10; %rad/s
omega3 = zeros(1, N);
bdot = zeros(1, N);
v0 = [0; 0]; %velocity of pin A
vB = zeros(2, N);
vP = zeros(2, N);

%variable for acceleration analysis
alpha2 = 0; %rad/s/s
alpha3 = zeros(1, N);
b_doubledot = zeros(1, N);
a0 = [0; 0]; %acceleration of crank pin A
aB = zeros(2, N);
aP = zeros(2, N);

%variable for force analysis
m2 = 0.090;
m3 = 0.240;
m4 = 0.200;
I2 = 0.0004;
I3 = 0.0003;
I4 = 0.0006;
xbar2 = [0.100; 0];
xbar3 = [0.100; 0];
xbar4 = [0.055; 0];

U2 = eye(2);
Z2 = zeros(2);
Z12 = zeros(1, 2);
Z21 = zeros(2, 1);
[FA, FB, FC, FD, FQ, FP] = deal(zeros(2, N));
M = zeros(1,N);


[v2, v3, v4] = deal(zeros(2, N));
[a2, a3, a4] = deal(zeros(2, N));
[P_Ext, P_Kin] = deal(zeros(1, N));

% pressure constants
D = 0.045; %piston diameter (m)
A = pi*D^2/4; %piston area
Patm = 101e3; %atmospheric pressure
V1 = A*(sqrt(a^2 + d^2) - p - H); %max cylinder volume
P = zeros(1, N); %cyclinder pressure 

% acceleration pulse parameters
theta2_min = 90*pi/180;   %crank position before compression
theta2_max = 20*pi/180;   %crank position at max compression
T = 2;  %time taken to fully compress
lambda = 2*pi/T;  %frequency
dt = T/(N-1); %time increment
t = 0:dt:T;
B = lambda*(theta2_max - theta2_min)/T;  %amplitude of accn pulse

% main loop
for i = 1:N
    theta2(i) = (B/lambda)*t(i) - (B/lambda^2)*sin(lambda*t(i)) + theta2_min;
    omega2 = (B/lambda)*(1 - cos(lambda*t(i)));
    alpha2 = B*sin(lambda*t(i));

    r = d - a * cos(theta2(i));
    s = a * sin(theta2(i));
    
    f = sqrt(r^2 + s^2);
    
    % open configuration - return positive value of b
    b(i) = c * cos(delta) + sqrt(f^2 - h^2);
    
    g = b(i) - c * cos(delta);
    
    theta3(i) = atan2((h*r - g*s), (g*r + h*s));
    
    theta4(i) = theta3(i) + delta;
    
    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));
    [e4, n4] = UnitVector(theta4(i));
    [eBP, nBP] = UnitVector(theta3(i) + gamma3);
    
    xB(:, i) = FindPos(x0, a, e2);
    xP(:, i) = FindPos(xB(:, i), p, e3);
    xC(:, i) = FindPos(xB(:, i), p, e3);
    xQ(:, i) = FindPos(x0, q, e2);

    % velocity analysis
    A_matrix = [(b(i)*n3 - c*n4), e3];
    b_vector = -a*omega2*n2;
    omega_vector = A_matrix\b_vector;
    omega3(i) = omega_vector(1);
    bdot(i) = omega_vector(2);

    %calculate point velocity
    vB(:, i) = FindVel(v0, a, omega2, n2);
    vP(:, i) = FindVel(vB(:, i), p, omega3(i), n3);

    % acceleration analysis
    C_matrix = [b(i)*n3-c*n4, e3];
    d_vector = -a*alpha2*n2 + a*omega2^2*e2 - 2*bdot(i)*omega3(i)*n3 + ...
    b(i)*omega3(i)^2*e3 - c*omega3(i)^2*e4;
    alpha_vector = C_matrix\d_vector;
    alpha3(i) = alpha_vector(1);
    b_doubledot(i) = alpha_vector(2);

    %calculate point acceleration
    aB(:, i) = FindAcc(a0, a, omega2, alpha2, e2, n2);
    aP(:, i) = FindAcc(aB(:, i), p, omega3(i), alpha3(i), eBP, nBP);

    % force analysis
    [eA2, nA2, LA2, s2A, s2B, s2Q] = LinkCG(a, q, 0, xbar2, theta2(i));
    [eB3, nB3, LB3, s3B, s3P, ~] = LinkCG(p, 0, gamma3, xbar3, theta3(i));
    [eD4, nD4, LD4, s4D, s4P, ~] = LinkCG(b(i) - p, 0, 0, xbar4, theta4(i));

    %rocker rotate at the same speed as slider
    omega4 = omega3(i);
    alpha4 = alpha3(i);

    v2(:, i) = FindVel(v0, LA2, omega2, nA2);
    v3(:, i) = FindVel(vB(:, i), LB3, omega3(i), nB3);
    v4(:, i) = FindVel(v0, LD4, omega4, nD4);
    vQ = FindVel(v0, q, omega2, n2);

    a2(:, i) = FindAcc(a0, LA2, omega2, alpha2, eA2, nA2);
    a3(:, i) = FindAcc(aB(:, i), LB3, omega3(i), alpha3(i), eB3, nB3);
    a4(:, i) = FindAcc(a0, LD4, omega4, alpha4, eD4, nD4);

    %pressure calculation
    V2 = A*(b(i) - p - H);
    P(i) = Patm*(V1/V2)^1.4;
    FP(:, i) = -A*(P(i) - Patm)*e3;

    S_matrix = [U2    -U2    Z2    Z2     U2    Z21;
                Z2     U2   -U2    Z2     Z2    Z21;
                Z2     Z2    U2   -U2     Z2    Z21;
                s2A'  -s2B'  Z12   Z12   s2Q'     0;
                Z12    s3B' -s3P'  Z12    Z12     1;
                Z12    Z12   s4P' -s4D'   Z12    -1;
                Z12    Z12    e3'  Z12    Z12     0;
                Z12    Z12    Z12  Z12    e2'     0];

    t_vector = [m2*a2(:, i);
                m3*a3(:, i) - FP(:, i);
                m4*a4(:, i);
                I2*alpha2;
                I3*alpha3(i) - dot(s3P, FP(:, i));
                I4*alpha4;
                0;
                0];
    
    f_vector = S_matrix\t_vector;

    FA(:, i) = [f_vector(1); f_vector(2)];
    FB(:, i) = [f_vector(3); f_vector(4)];
    FC(:, i) = [f_vector(5); f_vector(6)];
    FD(:, i) = [f_vector(7); f_vector(8)];
    FQ(:, i) = [f_vector(9); f_vector(10)];
    M(i) = f_vector(11);

    % Energy method validation
    P_FP = dot(FP(:, i), vP(:, i));
    P_FQ = dot(FQ(:, i), vQ);
    P_Ext(i) = P_FP + P_FQ;

    P2 = InertialPower(m2, I2, v2(:, i), a2(:, i), omega2, alpha2);
    P3 = InertialPower(m3, I3, v3(:, i), a3(:, i), omega3(i), alpha3(i));
    P4 = InertialPower(m4, I4, v4(:, i), a4(:, i), omega4, alpha4);
    P_Kin(i) = P2 + P3 + P4;
end

%--lines of code below are purely for validation of above point
%velocity. using finite difference method
timestep = dt;
vPx_estimated = FiniteDiffMethod(xP(1, :), timestep);
vPy_estimated = FiniteDiffMethod(xP(2, :), timestep);

%estimate alpha3 and d_doubledot with numerical method
alpha3_estimate = FiniteDiffMethod(omega3, timestep);
b_doubledot_estimate = FiniteDiffMethod(bdot, timestep);
aPx_estimate = FiniteDiffMethod(vP(1, :), timestep);
aPy_estimate = FiniteDiffMethod(vP(2, :), timestep);

%create a tiled layout
layout = tiledlayout(3, 2);

% plot pressure
nexttile(3);
Pgage = P-Patm;
plot(t, Pgage*1e-06);
set(gca, 'YTick', 0:0.2:2)
ylabel('Gage Pressure (Pa)');
legend('Pgage', Location='best');
grid on

% plot velocity or acceleration
nexttile(2);
plot(t, aP(1, :));
hold on
plot(t, aP(2, :));
%plot estimated value from finite difference
plot(t, aPy_estimate, '.');
ylabel('Velocity (m/s)');

% plot energy method validation
nexttile(4);
plot(t, P_Ext);
hold on
plot(t, P_Kin, '.');
ylabel('Power (W)');
legend('P\_Ext', 'P\_Kin', 'Location', 'best');

% plot pedal force
nexttile(5);
plot(t, FQ(1, :));
hold on
plot(t, FQ(2, :));
ylim([-500 200]);
ylabel('Pedal Force (N)');
legend('FQ\_x', 'FQ\_y', Location='best');
grid on

% plot mechanical advantage
nexttile(6);
plot(theta2*180/pi, vecnorm(FP) ./ vecnorm(FQ));
xlabel('Crank Angle (degree)');
ylabel('Mechanical Advantage, FP/FQ');
ylim([0 8])
grid on

% plot path
nexttile(1);
plot(xB(1, :), xB(2, :), '--');
hold on
%plot(xC(1, :), xC(2, :), '--');
plot(xP(1, :), xP(2, :), '--');
axis equal
ylim([0 0.25])

%initialise linkage plot
h1 = plot([x0(1), xB(1, 1), xQ(1, 1)], [x0(2), xB(2, 1), xQ(2, 1)]);
h3 = plot([xB(1), xC(1, 1), xP(1, 1)], [xB(2), xC(2, 1), xP(2, 1)]);
h4 = plot([xD(1), xP(1, 1)], [xD(2), xP(2, 1)]);

%initialise joint marker
h5 = plot([x0(1), xB(1, 1), xC(1, 1), xD(1), xP(1, 1), xQ(1, 1)], ...
    [x0(2), xB(2, 1), xC(2, 1), xD(2), xP(2, 1), xQ(2, 1)], 'o', 'MarkerFaceColor', 'k');

%initialise label
t1 = text(x0(1) + 0.005, x0(2), 0, "A");
t2 = text(xB(1, 1) + 0.005, xB(2, 1), 0, "B");
%t3 = text(xC(1, 1) + 0.005, xC(2, 1), 0, "C");
t4 = text(xP(1, 1) + 0.005, xP(2, 1), 0, "P");
t5 = text(xD(1) + 0.005, xD(2), 0, "D");
t6 = text(xQ(1, 1) + 0.005, xQ(2, 1), 0, "Q");

% animation loop
j = 1;
direction = 1;
while j <= N
    if ~isgraphics(layout)
        break
    end
    set(h1, 'XData', [x0(1), xB(1, j), xQ(1, j)], 'YData', [x0(2), xB(2, j), xQ(2, j)]);
    set(h3, 'XData', [xB(1, j), xC(1, j), xP(1, j)], 'YData', [xB(2, j), xC(2, j), xP(2, j)]);
    set(h4, 'XData', [xD(1), xP(1, j)], 'YData', [xD(2), xP(2, j)]);

    %update joint marker
    set(h5, 'XData', [x0(1), xB(1, j), xC(1, j), xD(1), xP(1, j), xQ(1, j)], ...
    'YData', [x0(2), xB(2, j), xC(2, j), xD(2), xP(2, j), xQ(2, j)])

    drawnow

    %update label position
    %check if text object still exist first, before update
    if all(isgraphics([t2, t4]))
        set(t2, 'Position', [xB(1, j) + 0.005, xB(2, j)]);
        %set(t3, 'Position', [xC(1, j) + 0.005, xC(2, j)]);
        set(t4, 'Position', [xP(1, j) + 0.005, xP(2, j)]);
        set(t6, 'Position', [xQ(1, j) + 0.005, xQ(2, j)]);
    end

    % if not grashoff, check if each extreme position has been reached. 
    % if yes, then change direction. 
    % if grashoff, use the logic that reset j = 1 every complete loop
    if grashoff == false
        if (j >= N && direction > 0) || (j <= 1 && direction < 0)
            direction = direction * -1;
        end
    elseif grashoff == true
        if j == N & loop == true
            j = 1;
        end
    end
    j = j + 1*direction;


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