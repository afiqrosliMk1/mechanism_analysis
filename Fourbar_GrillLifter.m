% Saved from Fourbar_Force_Analysis.m
% Afiq Rosli, 19/10/2025
% the methods of finding coupler angle (theta 3) and rocker angle (theta 4)
% is taken from Karl B Dyer et al

% prepare workspace
clear variables; close all; clc;
% crank length (m)
a = 0.03004;
% coupler length (m)
b = 0.06000;
% rocker length (m)
c = 0.03870;
% ground link length (m)
d = 0.01481;
% length from B to P (m)
p = 0.59730; 
% angle between BP and coupler (radian)
gamma = -55 *pi/180; 

grashoff = false;
% check for grashoff condition
[S, L, P, Q] = slpq(a , b, c, d);
if S + L < P + Q
    grashoff = true;
    fprintf('grashoff\n');
else
    fprintf('non-grashoff\n');
end

% initialise variables
loop = true;
N = 361;
theta2 = zeros(1, N);
theta3 = zeros(1, N);
theta3_prime = zeros(1, N);
theta4 = zeros(1, N);
delta = 0;
r = 0;
s = 0;
f = 0;
g = 0;
h = 0;
x0 = [0; 0];

theta2_min = 87.04*pi/180;
theta2_max = 179.78*pi/180;

xD = [d; 0];
xB = zeros(2, N);
xC = zeros(2, N);
xP = zeros(2, N);

%store velocity analysis result
%omega2 = 10; %rad/s
omega2 = zeros(1, N);
omega3 = zeros(1, N);
omega4 = zeros(1, N);

v0 = [0; 0]; % velocity at crank pivot
vB = zeros(2, N);
vP = zeros(2, N);

[v2, v3, v4] = deal(zeros(2, N));

% acceleration analysis variables
%alpha2 = 0;
alpha2 = zeros(1, N);
alpha3 = zeros(1, N);
alpha4 = zeros(1, N);

a0 = [0; 0]; % acceleration at crank pivot
aB = zeros(2, N);
aP = zeros(2, N);
aC = zeros(2, N);

%force analysis variables
m2 = 0.0; % crank mass [kg]
m3 = 60; % coupler mass [kg]
m4 = 0.0; % rocker mass [kg]
I2 = 0.0; % crank moment of inertia [kg.m^2]
I3 = 1.35; % coupler moment of inertia [kg.m^2]
I4 = 0.0; % rocker momentof inertia [kg.m^2]
gravity = 9.81; % gravitational acceleration [m/s^2]
Fg = [-m3*gravity; 0]; % force of gravity on the lid
xk0 = 0.01; % unstretched length of the spring 
k = 90000; % spring constant [N/m]

xbar2 = [a/2; 0];
xbar3 = [0.2878; -0.1745];
xbar4 = [c/2; 0];

U2 = eye(2);
Z2 = zeros(2);
Z21 = zeros(2, 1);
Z12 = zeros(1, 2);

[a2, a3, a4] = deal(zeros(2, N)); % accn of COG of each link [m/s^2]

[FA, FB, FC, FD, FP] = deal(zeros(2, N));
T3 = zeros(1, N);

[P_Ext, P_Kin] = deal(zeros(1, N));

T = 2; % secs
dt = T/(N-1);
t = 0:dt:T;
lambda = 2*pi/T;
A = lambda*(theta2_max - theta2_min)/T;

% main loop
for i = 1:N

    theta2(i) = (A/lambda)*t(i) - (A/lambda^2)*sin(lambda*t(i)) + theta2_min;
    omega2(i) = (A/lambda)*(1 - cos(lambda*t(i)));
    alpha2(i) = A*sin(lambda*t(i));

    r = d - a * cos(theta2(i));
    s = a * sin(theta2(i));

    f = sqrt(r^2 + s^2);
    
    % open configuration
    cos_val = (b^2 + c^2 - f^2) / (2*b*c);
    
    % Clamp to [-1, 1]. Important to handle floating point error. otherwise
    % delta will give you 0 + 2.107342425544702e-08i for theta2_min angle.
    cos_val = min(1, max(-1, cos_val));   
    delta = acos(cos_val);
    
    g = b - c * cos(delta);
    h = c * sin(delta);

    theta3(i) = atan2((h*r - g*s), (g*r + h*s));
    theta4(i) = theta3(i) + delta;

    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));
    [e4, n4] = UnitVector(theta4(i));

    [eBP, nBP] = UnitVector(theta3(i) + gamma);

    xB(:, i) = FindPos(x0, a, e2);
    xC(:, i) = FindPos(xD, c, e4);
    xP(:, i) = FindPos(xB(:, i), p, eBP);

    % velocity analysis
    A_matrix = [b*n3, -c*n4];
    b_vector = -a*omega2(i)*n2;
    omega_vector = A_matrix \ b_vector;
    omega3(i) = omega_vector(1);
    omega4(i) = omega_vector(2);

    %solve for point velocity
    vB(:, i) = FindVel(v0, a, omega2(i), n2);
    vP(:, i) = FindVel(vB(:, i), p, omega3(i), nBP);

    % acceleration analysis
    C_matrix = [b*n3, -c*n4];
    d_vector = -alpha2(i)*a*n2 + omega2(i)^2*a*e2 + omega3(i)^2*b*e3 - omega4(i)^2*c*e4;
    %dr constans somehow use the one below, omitting the -alpha2*a*n2 term
    %d_vector = omega2(i)^2*a*e2 + omega3(i)^2*b*e3 - omega4(i)^2*c*e4;
    %pending his confirmation. email sent 23/10/2025
    %dr constans replied on 26/10/2025 confirming that omitting
    %-alpha2*a*n2 was an error
    alpha_matrix = C_matrix\d_vector;
    alpha3(i) = alpha_matrix(1);
    alpha4(i) = alpha_matrix(2);

    %solve for point acceleration
    aB(:, i) = FindAcc(a0, a, omega2(i), alpha2(i), e2, n2);
    aP(:, i) = FindAcc(aB(:, i), p, omega3(i), alpha3(i), eBP, nBP);
    aC(:, i) = FindAcc(a0, c, omega4(i), alpha4(i), e4, n4);

    %force analysis
    [eA2, nA2, LA2, s2A, s2B, ~] = LinkCG(a, 0, 0, xbar2, theta2(i));
    [eB3, nB3, LB3, s3B, s3C, s3P] = LinkCG(b, p, gamma, xbar3, theta3(i));
    [eD4, nD4, LD4, s4D, s4C, ~] = LinkCG(c, 0, 0, xbar4, theta4(i));

    a2(:, i) = FindAcc(a0, LA2, omega2(i), alpha2(i), eA2, nA2);
    a3(:, i) = FindAcc(aB(:, i), LB3, omega3(i), alpha3(i), eB3, nB3);
    a4(:, i) = FindAcc(a0, LD4, omega4(i), alpha4(i), eD4, nD4);

    xQ = [-0.05; 0];
    rBQ = xQ - xB(:, i);
    eBQ = (1/norm(rBQ))*rBQ;

    rCP = xP(:, i) - xC(:, i);
    eCP = (1/norm(rCP))*rCP;

    Fk = k*(norm(xB(:, i) - xQ) - xk0) * eBQ;

    S_matrix = [U2    -U2    Z2    Z2    Z2;
                Z2     U2   -U2    Z2    U2;
                Z2     Z2    U2   -U2    Z2;
                s2A'  -s2B'  Z12   Z12   Z12;
                Z12    s3B' -s3C'  Z12   s3P';
                Z12    Z12   s4C' -s4D'  Z12;
                Z12    Z12   Z12   Z12   eCP'];

    t_vector = [m2*a2(:, i);
                m3*a3(:, i) - Fk - Fg;
                m4*a4(:, i);
                I2*alpha2(i);
                I3*alpha3(i) - dot(Fk, s3B);
                I4*alpha4(i);
                0];

    f_vector = S_matrix\t_vector;

    FA(:, i) = [f_vector(1); f_vector(2)];
    FB(:, i) = [f_vector(3); f_vector(4)];
    FC(:, i) = [f_vector(5); f_vector(6)];
    FD(:, i) = [f_vector(7); f_vector(8)];
    FP(:, i) = [f_vector(9); f_vector(10)];

    T3(i) = dot(FP(:, i), s3P);

    % Energy method for validation
    v2(:, i) = FindVel(v0, LA2, omega2(i), nA2);
    v3(:, i) = FindVel(vB(:, i), LB3, omega3(i), nB3);
    v4(:, i) = FindVel(v0, LD4, omega4(i), nD4);

    P_FP = dot(FP(:, i), vP(:, i));
    P_Fg = dot(Fg, v3(:, i));
    P_Fk = dot(Fk, vB(:, i));
    P_Ext(i) = P_FP + P_Fg + P_Fk;

    P2 = InertialPower(m2, I2, v2(:, i), a2(:, i), omega2(i), alpha2(i));
    P3 = InertialPower(m3, I3, v3(:, i), a3(:, i), omega3(i), alpha3(i));
    P4 = InertialPower(m4, I4, v4(:, i), a4(:, i), omega4(i), alpha4(i));
    P_Kin(i) = P2 + P3 + P4;
end

%%%-- finite diff method. purely for validation of velocity result-%%
%timestep = 2*pi/((N-1)*omega2);
timestep = dt;
vPx_estimated = FiniteDiffMethod(xP(1, :), timestep);
vPy_estimated = FiniteDiffMethod(xP(2, :), timestep);
omega3_estimated = FiniteDiffMethod(theta3(1, :), timestep);
omega4_estimated = FiniteDiffMethod(theta4(1, :), timestep);

% estimate acceleration with numerical method
alpha3_estimate = FiniteDiffMethod(omega3, timestep);
alpha4_estimate = FiniteDiffMethod(omega4, timestep);
aPx_estimate = FiniteDiffMethod(vP(1, :), timestep);
aPy_estimate = FiniteDiffMethod(vP(2, :), timestep);

%%PLOTTING RESULT ON A TILED LAYOUT
layout = tiledlayout(3, 2);

% plot location
nexttile(2);
plot(t, P_Ext);
hold on;
plot(t, P_Kin, '.');
legend('P\_Ext', 'P\_Kin', Location='best');
ylabel('Power');
xlabel('time (s)');

% plot velocity
nexttile(3);
plot(t, vP(1, :));
hold on
plot(t, vP(2, :));
% --plot below purely to validate the analytical result vs estimated
plot(t, vPy_estimated, '.');
legend('vPx\_analytical', 'vPy\_analytical', 'vPy\_estimated', 'Location','best');
ylabel('Velocity');
xlabel('time (s)');

% plot acceleration
nexttile(5);
plot(t, aP(2, :));
hold on;
plot(t, aP(1, :));
plot(t, aPy_estimate, '.');
xlabel('time (s)');
ylabel('Acceleration');
legend('aPy', 'aPy', 'aPy\_estimated', 'Location', 'best');

% plot torque
nexttile(4);
plot(t, T3);
xlabel('time (s)');
ylabel('Opening Torque (Nm)');
grid on;

% plot power
nexttile(6);
plot(t, FP(1, :));
hold on
plot(t, FP(2, :));
ylabel('User Force, FP (N)');
xlabel('time (s)');
legend('FP\_x', 'FP\_y', 'Location', 'best');

%plot non moving linkage
nexttile(1);
plot(xB(1, :), xB(2, :));
hold on
plot(xC(1, :), xC(2, :));
plot(xP(1, :), xP(2, :));
axis equal

plot([x0(1), xB(1, 1), xC(1, 1)], [x0(2), xB(2, 1), xC(2, 1)], HandleVisibility='off');
plot([xB(1, 1), xP(1, 1), xC(1, 1)], [xB(2, 1), xP(2, 1), xC(2, 1)], HandleVisibility='off');
plot([xD(1), xC(1, 1)], [xD(2), xC(2, 1)], HandleVisibility='off');

%initialise linkage plot
h1 = plot([x0(1), xB(1, 1), xC(1, 1)], [x0(2), xB(2, 1), xC(2, 1)], HandleVisibility='off');
h2 = plot([xB(1, 1), xP(1, 1), xC(1, 1)], [xB(2, 1), xP(2, 1), xC(2, 1)], HandleVisibility='off');
h3 = plot([xD(1), xC(1, 1)], [xD(2), xC(2, 1)], HandleVisibility='off');

%h4 = patch([xB(1, 1), xP(1, 1), xC(1, 1)], [xB(2, 1), xP(2, 1), xC(2, 1)], [0.5, 0.5, 1], 'FaceAlpha', 0.3);

%initialise joint marker
h5 = plot([x0(1), xB(1, 1), xC(1, 1), xD(1), xP(1, 1)], ...
    [x0(2), xB(2, 1), xC(2, 1), xD(2), xP(2, 1)], 'o', 'MarkerFaceColor', 'k');

%initialise label
t1 = text(x0(1) + 0.005, x0(2), 0, "A");
t2 = text(xB(1, 1) + 0.005, xB(2, 1), 0, "B");
t3 = text(xC(1, 1) + 0.005, xC(2, 1), 0, "C");
t4 = text(xP(1, 1) + 0.005, xP(2, 1), 0, "P");
t5 = text(xD(1) + 0.005, xD(2), 0, "D");

% animation loop
j = 1;
direction = 1;
while j <= N
    if ~ishandle(layout)
        break
    end

    %update linkage
    set(h1, 'XData', [x0(1), xB(1, j), xC(1, j)], 'YData', [x0(2), xB(2, j), xC(2, j)]);
    set(h2, 'XData', [xB(1, j), xP(1, j), xC(1, j)], 'YData', [xB(2, j), xP(2, j), xC(2, j)]);
    set(h3, 'XData', [xD(1), xC(1, j)], 'YData', [xD(2), xC(2, j)]);

    %update joint marker
    set(h5, 'XData', [x0(1), xB(1, j), xC(1, j), xD(1), xP(1, j)], ...
    'YData', [x0(2), xB(2, j), xC(2, j), xD(2), xP(2, j)])

    %update coupler patch
    %set(h4, 'XData', [xB(1, j), xP(1, j), xC(1, j)], 'YData', [xB(2, j), xP(2, j), xC(2, j)]);
    drawnow

    %update label position
    %check if text object still exist first, before update
    if all(isgraphics([t2, t3, t4]))
        set(t2, 'Position', [xB(1, j) + 0.005, xB(2, j)]);
        set(t3, 'Position', [xC(1, j) + 0.005, xC(2, j)]);
        set(t4, 'Position', [xP(1, j) + 0.005, xP(2, j)]);
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