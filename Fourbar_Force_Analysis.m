% Saved from Fourbar_Position_Analysis_with_non_grashoff.m
% ++with solver for non grashoff case
% conducts a force analysis of the fourbar linkage
% positions of points B, C and P
% Afiq Rosli, 19/10/2025
% the methods of finding coupler angle (theta 3) and rocker angle (theta 4)
% is taken from Karl B Dyer et al

% prepare workspace
clear variables; close all; clc;
% crank length (m)
a = 0.13;
% coupler length (m)
b = 0.2;
% rocker length (m)
c = 0.17;
% ground link length (m)
d = 0.22;
% length from B to P (m)
p = 0.150; 
% angle between BP and coupler (radian)
gamma = 20 * pi / 180; 

% if grashoff condition met, grashoff = true, otherwise false
grashoff = true;
% linkage class. initialise to 0.
class = 0;
% driving link clocwise or counterclockwise
ccw = true;

% check for grashoff condition
[S, L, P, Q] = slpq(a , b, c, d);
if S + L < P + Q
    if S == d
        % short link is ground
        class = 1;
    elseif S == a
        % short link is crank
        class = 2;
    elseif S == b
        % short link is coupler
        class = 3;
    elseif S == c
        % short link is rocker
        class = 4;
    end
    fprintf("grashoff class %d\n", class)
    theta2_min = 0;
    theta2_max = 2*pi;
else
    if L == a
        fprintf("non grashoff class 7 - longest link is crank (link ""a"")\n")
        class = 7;
    elseif L == b
        fprintf("non grashoff class 8&9 - longest link is coupler (link ""b"")\n")
        class = 8;
        cw = false;
        theta2_max = acos((a^2 + d^2 - (b - c)^2) / (2*a*d));
        theta2_min = -theta2_max;
    elseif L == c
        fprintf("non grashoff class 10 - longest link is rocker (link ""c"")\n")
        class = 10;
    elseif L == d
        fprintf("non grashoff class 5&6 - longest link is ground (link ""d"")\n")
        class = 5;
        theta2_max = acos((a^2 + d^2 - (b + c)^2) / (2*a*d));
        theta2_min = -theta2_max;
    end
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
if class == 2 || class == 1
    xD = [d; 0];
elseif class == 4
    xD = [-d; 0];
else
    fprintf("warning, xD not defined\n")
end
xB = zeros(2, N);
xC = zeros(2, N);
xP = zeros(2, N);

%store velocity analysis result
omega2 = 10; %rad/s
omega3 = zeros(1, N);
omega4 = zeros(1, N);

v0 = [0; 0]; % velocity at crank pivot
vB = zeros(2, N);
vP = zeros(2, N);

[v2, v3, v4] = deal(zeros(2, N));

% acceleration analysis variables
alpha2 = 0;
alpha3 = zeros(1, N);
alpha4 = zeros(1, N);

a0 = [0; 0]; % acceleration at crank pivot
aB = zeros(2, N);
aP = zeros(2, N);
aC = zeros(2, N);

%force analysis variables
m2 = 0.124; % crank mass [kg]
m3 = 0.331; % coupler mass [kg]
m4 = 0.157; % rocker mass [kg]
I2 = 0.000255; % crank moment of inertia [kg.m^2]
I3 = 0.001188; % coupler moment of inertia [kg.m^2]
I4 = 0.000503; % rocker momentof inertia [kg.m^2]

xbar2 = [a/2; 0];
xbar3 = [0.1071; 0.0148];
xbar4 = [c/2; 0];

U2 = eye(2);
Z2 = zeros(2);
Z21 = zeros(2, 1);
Z12 = zeros(1, 2);

[a2, a3, a4] = deal(zeros(2, N)); % accn of COG of each link [m/s^2]

[FA, FB, FC, FD] = deal(zeros(2, N));
T2 = zeros(1, N);

FP = [0; -100]; % force at point P [N]
T4 = 0; % torque applied to rocker [N.m]

[P_Ext, P_Kin] = deal(zeros(1, N));

% main loop
for i = 1:N
    if ccw == true
        % if driving link rotates ccw, then theta2 takes positive angle
        theta2(i) = (i - 1) * (theta2_max - theta2_min) / (N - 1) + theta2_min;
    elseif ccw == false
        % if driving link rotates cw, then theta2 takes negative angle
        theta2(i) = -1 * (i - 1) * (2*pi - (theta2_max - theta2_min)) / (N - 1) + theta2_min;
    end
  
    if class == 2 || class == 1
        r = d - a * cos(theta2(i));
        s = a * sin(theta2(i));
    elseif class == 4
        r = d + c * cos(theta2(i));
        s = c * sin(theta2(i));
    else
        fprintf("class other than 2 and 4 not coverted.. return")
        return
    end

    f = sqrt(r^2 + s^2);
    
    % open configuration
    if class == 2 || class == 1
        cos_val = (b^2 + c^2 - f^2) / (2*b*c);
    elseif class == 4
        cos_val = (a^2 + b^2 - f^2) / (2*a*b);
    end
    
    % Clamp to [-1, 1]. Important to handle floating point error. otherwise
    % delta will give you 0 + 2.107342425544702e-08i for theta2_min angle.
    cos_val = min(1, max(-1, cos_val));   
    delta = acos(cos_val);
    
    if class == 2 || class == 1
        g = b - c * cos(delta);
        h = c * sin(delta);
    elseif class == 4
        g = b - a * cos(delta);
        h = a * sin(delta);
    end

    if class == 2 || class == 1
        theta3(i) = atan2((h*r - g*s), (g*r + h*s));
        theta4(i) = theta3(i) + delta;
    elseif class == 4
        % class 4 is the mirror of class 2
        theta3_prime(i) = atan2((h*r - g*s), (g*r + h*s));
        theta3(i) = pi - theta3_prime(i);
        theta4(i) = theta3(i) - delta;
    end

    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));
    [e4, n4] = UnitVector(theta4(i));

    if class == 2 || class == 1
        [eBP, nBP] = UnitVector(theta3(i) + gamma);
    elseif class == 4
        [eCP, nCP] = UnitVector(theta3(i) - gamma);
    end

    if class == 2 || class == 1
        xB(:, i) = FindPos(x0, a, e2);
        xC(:, i) = FindPos(xD, c, e4);
        xP(:, i) = FindPos(xB(:, i), p, eBP);
    elseif class == 4
        xB(:, i) = FindPos(xD, a, e4);
        xC(:, i) = FindPos(x0, c, e2);
        xP(:, i) = FindPos(xC(:, i), p, eCP);
    end

    % velocity analysis
    A_matrix = [b*n3, -c*n4];
    b_vector = -a*omega2*n2;
    omega_vector = A_matrix \ b_vector;
    omega3(i) = omega_vector(1);
    omega4(i) = omega_vector(2);

    %solve for point velocity
    vB(:, i) = FindVel(v0, a, omega2, n2);
    vP(:, i) = FindVel(vB(:, i), p, omega3(i), nBP);

    % acceleration analysis
    C_matrix = [b*n3, -c*n4];
    d_vector = -alpha2*a*n2 + omega2^2*a*e2 + omega3(i)^2*b*e3 - omega4(i)^2*c*e4;
    alpha_matrix = C_matrix\d_vector;
    alpha3(i) = alpha_matrix(1);
    alpha4(i) = alpha_matrix(2);

    %solve for point acceleration
    aB(:, i) = FindAcc(a0, a, omega2, alpha2, e2, n2);
    aP(:, i) = FindAcc(aB(:, i), p, omega3(i), alpha3(i), eBP, nBP);
    aC(:, i) = FindAcc(a0, c, omega4(i), alpha4(i), e4, n4);

    %force analysis
    [eA2, nA2, LA2, s2A, s2B, ~] = LinkCG(a, 0, 0, xbar2, theta2(i));
    [eB3, nB3, LB3, s3B, s3C, s3P] = LinkCG(b, p, gamma, xbar3, theta3(i));
    [eD4, nD4, LD4, s4D, s4C, ~] = LinkCG(c, 0, 0, xbar4, theta4(i));

    a2(:, i) = FindAcc(a0, LA2, omega2, alpha2, eA2, nA2);
    a3(:, i) = FindAcc(aB(:, i), LB3, omega3(i), alpha3(i), eB3, nB3);
    a4(:, i) = FindAcc(a0, LD4, omega4(i), alpha4(i), eD4, nD4);


    S_matrix = [U2    -U2    Z2    Z2    Z21;
                Z2     U2   -U2    Z2    Z21;
                Z2     Z2    U2   -U2    Z21;
                s2A'  -s2B'  Z12   Z12   1;
                Z12    s3B' -s3C'  Z12   0;
                Z12    Z12   s4C' -s4D'  0];

    t_vector = [m2*a2(:, i);
                m3*a3(:, i) - FP;
                m4*a4(:, i);
                I2*alpha2;
                I3*alpha3(i) - dot(s3P, FP);
                I4*alpha4(i) - T4];

    f_vector = S_matrix\t_vector;

    FA(:, i) = [f_vector(1); f_vector(2)];
    FB(:, i) = [f_vector(3); f_vector(4)];
    FC(:, i) = [f_vector(5); f_vector(6)];
    FD(:, i) = [f_vector(7); f_vector(8)];
    T2(i) = f_vector(9);

    % Energy method for validation
    v2(:, i) = FindVel(v0, LA2, omega2, nA2);
    v3(:, i) = FindVel(vB(:, i), LB3, omega3(i), nB3);
    v4(:, i) = FindVel(v0, LD4, omega4(i), nD4);

    P_FP = dot(FP, vP(:, i));
    P_T2 = T2(i)*omega2;
    P_T4 = T4*omega4(i);
    P_Ext(i) = P_T2 + P_FP + P_T4;

    P2 = InertialPower(m2, I2, v2(:, i), a2(:, i), omega2, alpha2);
    P3 = InertialPower(m3, I3, v3(:, i), a3(:, i), omega3(i), alpha3(i));
    P4 = InertialPower(m4, I4, v4(:, i), a4(:, i), omega4(i), alpha4(i));
    P_Kin(i) = P2 + P3 + P4;
end

%%%-- finite diff method. purely for validation of velocity result-%%
timestep = 2*pi/((N-1)*omega2);
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
plot(theta2*180/pi, xP(1, :));
xlim([0 360]);
set(gca, 'xtick', 0:60:360);
xlabel('crank angle');
ylabel('Position');
legend('xP', Location='best');

% plot velocity
nexttile(3);
plot(theta2*180/pi, vP(2, :));
xlim([0, 360]);
set(gca, 'xtick', 0:60:360)

% --plot below purely to validate the analytical result vs estimated
% -- change the plot to omega3, omega4, alpha3, alpha4 or vP as you wish
hold on
plot(theta2*180/pi, vPy_estimated, '.');
legend('analytical', 'estimated', 'Location','best');
xlabel('crank angle');
ylabel('Velocity');

% plot acceleration
nexttile(5);
plot(theta2*180/pi, aP(2, :));
xlim([0, 360]);
set(gca, 'xtick', 0:60:360)
xlabel('crank angle');
ylabel('Acceleration');
legend('aP', Location='best');

% plot torque
nexttile(4);
plot(theta2*180/pi, T2);
xlim([0 360]);
set(gca, 'xtick', 0:60:360);
xlabel('crank angle');
ylabel('driving torque [N.m]');
legend('T2', Location='best');
grid on;

% plot power
nexttile(6);
plot(theta2*180/pi, P_Ext);
hold on;
plot(theta2*180/pi, P_Kin, '.');
xlim([0 360]);
legend('P\_Ext', 'P\_Kin', Location='best');
xlabel('crank angle');
ylabel('Power');

% plot path
nexttile(1);
plot(xB(1, :), xB(2, :));
hold on
plot(xC(1, :), xC(2, :));
plot(xP(1, :), xP(2, :));
axis equal

%initialise linkage plot
if class == 2 || class == 1
    h1 = plot([x0(1), xB(1, 1), xC(1, 1)], [x0(2), xB(2, 1), xC(2, 1)], HandleVisibility='off');
    h2 = plot([xB(1, 1), xP(1, 1), xC(1, 1)], [xB(2, 1), xP(2, 1), xC(2, 1)], HandleVisibility='off');
    h3 = plot([xD(1), xC(1, 1)], [xD(2), xC(2, 1)], HandleVisibility='off');
elseif class == 4
    h1 = plot([x0(1), xC(1, 1), xB(1, 1)], [x0(2), xC(2, 1), xB(2, 1)]);
    h2 = plot([xC(1, 1), xP(1, 1), xB(1, 1)], [xC(2, 1), xP(2, 1), xB(2, 1)]);
    h3 = plot([xD(1), xB(1, 1)], [xD(2), xB(2, 1)]);
end

h4 = patch([xB(1, 1), xP(1, 1), xC(1, 1)], [xB(2, 1), xP(2, 1), xC(2, 1)], [0.5, 0.5, 1]);

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

    if class == 2 || class == 1
        set(h1, 'XData', [x0(1), xB(1, j), xC(1, j)], 'YData', [x0(2), xB(2, j), xC(2, j)]);
        set(h2, 'XData', [xB(1, j), xP(1, j), xC(1, j)], 'YData', [xB(2, j), xP(2, j), xC(2, j)]);
        set(h3, 'XData', [xD(1), xC(1, j)], 'YData', [xD(2), xC(2, j)]);
    elseif class == 4
        set(h1, 'XData', [x0(1), xC(1, j), xB(1, j)], 'YData', [x0(2), xC(2, j), xB(2, j)]);
        set(h2, 'XData', [xC(1, j), xP(1, j), xB(1, j)], 'YData', [xC(2, j), xP(2, j), xB(2, j)]);
        set(h3, 'XData', [xD(1), xB(1, j)], 'YData', [xD(2), xB(2, j)]);
    end

    %update joint marker
    set(h5, 'XData', [x0(1), xB(1, j), xC(1, j), xD(1), xP(1, j)], ...
    'YData', [x0(2), xB(2, j), xC(2, j), xD(2), xP(2, j)])

    %update coupler patch
    set(h4, 'XData', [xB(1, j), xP(1, j), xC(1, j)], 'YData', [xB(2, j), xP(2, j), xC(2, j)]);
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
    if class == 5
        if (j >= N && direction > 0) || (j <= 1 && direction < 0)
            direction = direction * -1;
        end
    elseif class == 6
        fprintf("under construction\n")
    elseif class == 7
        if (j >= N && direction > 0) || (j <= 1 && direction < 0)
            direction = direction * -1;
        end
    elseif class == 8
        fprintf("under construction\n")
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