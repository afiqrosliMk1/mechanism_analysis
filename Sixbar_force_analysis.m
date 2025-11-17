% sixbar_Force_Analysis.m
% conducts a force analysis on the Stephenson Type I linkage
% Afiq Rosli (16/11/2025)

% prepare workspace
clear variables; close all; clc;
a = 0.070; % crank length (m)
b = 0.100; % coupler length (m)
c = 0.09;  % rocker length (m)
d = 0.110; % ground link length (m)
p = 0.150; % length from A to E (m)
q = 0.150; % length from D to F (m)
u = 0.120; % length from E to G (m)
v = 0.160; % length from G to F (m)
w = 0.150;
gamma2 = 20*pi/180; % angle between AE and a (radian)
gamma4 = -20*pi/180;
gamma6 = -20*pi/180;

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
[xB, xC, xE, xF, xG] = deal(zeros(2, N));
alpha = 0;
beta = 0;
theta2_prime = 0;
theta5_prime = 0;
theta6_prime = 0;
theta5 = zeros(1, N);
theta6 = zeros(1, N);

% for velocity analysis
omega2 = 10; %rad/s
omega_vector = zeros(1, N);
[omega3, omega4, omega5, omega6] = deal(zeros(1, N));
v0 = [0; 0]; %velocity of crank pin
[vB, vE, vF, vG, vW] = deal(zeros(2, N));
[v2, v3, v4, v5, v6] = deal(zeros(2, N));

% for estimated velocity
[omega3_estimated, omega4_estimated, omega5_estimated, omega6_estimated] = deal(zeros(1, N));

% variables for acceleration analysis
alpha2 = 0; %rad/s/s
[alpha3, alpha4, alpha5, alpha6] = deal(zeros(1, N));
a0 = [0; 0]; %accelaration of crank pin A
[aB, aE, aF, aG] = deal(zeros(2, N));
[a2, a3, a4, a5, a6] = deal(zeros(2, N));

% variables for force analysis
U2 = eye(2); Z2 = zeros(2); Z21 = zeros(2, 1); Z12 = zeros(1, 2);

m2 = 0.200; m3 = 0.100; m4 = 0.200; m5 = 0.100; m6 = 0.150;
I2 = 0.0002; I3 = 0.0001; I4 = 0.0002; I5 = 0.0001; I6 = 0.00015;
xbar2 = [0.035; 0.0175];
xbar3 = [0.050; 0.0];
xbar4 = [0.045; -0.0225];
xbar5 = [0.060; 0.0];
xbar6 = [0.080; 0.04];

FW = [0; -100];
[FA, FB, FC, FD, FE, FF, FG] = deal(zeros(2, N));
T2 = zeros(1, N);
[P_Kin, P_Ext] = deal(zeros(1, N));

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
        [eAE, nAE] = UnitVector(theta2(i) + gamma2);
        [eDF, nDF] = UnitVector(theta4(i) + gamma4);
    elseif class == 4
        [eCP, nCP] = UnitVector(theta2(i) - gamma2);
    end

    if class == 2 || class == 1
        xB(:, i) = FindPos(x0, a, e2);
        xC(:, i) = FindPos(xD, c, e4);
        xE(:, i) = FindPos(x0, p, eAE);
        xF(:, i) = FindPos(xD, q, eDF);
    elseif class == 4
        xB(:, i) = FindPos(xD, a, e4);
        xC(:, i) = FindPos(x0, c, e2);
        xE(:, i) = FindPos(xC(:, i), p, eCP);
    end

    % find alpha and beta - sixbar linkage
    alpha = atan2((xE(2, i) - xB( 2, i)), (xE(1, i) - xB(1, i)));
    beta = atan2((xF(2, i) - xB(2 , i)), (xF(1, i) - xB(1, i)));

    a_prime = sqrt((xE(2, i) - xB(2, i))^2 + (xE(1, i) - xB(1, i))^2);
    d_prime = sqrt((xF(2, i) - xB(2, i))^2 + (xF(1, i) - xB(1, i))^2);

    theta2_prime = alpha - beta;
    r = d_prime - a_prime * cos(theta2_prime);
    s = a_prime * sin(theta2_prime);

    f = sqrt(r^2 + s^2);

    cos_val = (u^2 + v^2 - f^2) / (2*u*v);
    cos_val = min(1, max(-1, cos_val));   
    delta = acos(cos_val);

    g = u - v * cos(delta);
    h = v * sin(delta);

    theta5_prime = atan2((h*r - g*s), (g*r + h*s));
    theta6_prime = theta5_prime + delta;

    theta5(i) = theta5_prime + beta;
    theta6(i) = theta6_prime + beta;

    [eEG, nEG] = UnitVector(theta5(i));
    xG(:, i) = FindPos(xE(:, i), u, eEG);

    % velocity analysis
    [e5, n5] = UnitVector(theta5(i));
    [e6, n6] = UnitVector(theta6(i));
    [eFW, nFW] = UnitVector(theta6(i) + gamma6);
 
    A_matrix = [b*n3, -c*n4, Z21, Z21;
                Z21, -q*nDF, u*n5, -v*n6];
    b_vector = [-a*omega2*n2; -p*omega2*nAE];
    omega_vector = A_matrix \ b_vector;
    omega3(i) = omega_vector(1);
    omega4(i) = omega_vector(2);
    omega5(i) = omega_vector(3);
    omega6(i) = omega_vector(4);

    vB(:, i) = FindVel(v0, a, omega2, n2);
    vE(:, i) = FindVel(v0, p, omega2, nAE);
    vF(:, i) = FindVel(v0, q, omega4(i), nDF);
    vG(:, i) = FindVel(vE(:, i), u, omega5(i), n5);
    vW(:, i) = FindVel(vF(:, i), w, omega6(i), nFW);

    % acceleration analysis
    C_matrix = [b*n3, -c*n4, Z21, Z21;
                Z21, -q*nDF, u*n5, -v*n6];
    d_vector = [-a*alpha2*n2 + a*omega2^2*e2 + b*omega3(i)^2*e3 - c*omega4(i)^2*e4;
                -p*alpha2*nAE + p*omega2^2*eAE + u*omega5(i)^2*e5 - v*omega6(i)^2*e6 - q*omega4(i)^2*eDF];
    alpha_vector = C_matrix\d_vector;
    alpha3(i) = alpha_vector(1);
    alpha4(i) = alpha_vector(2);
    alpha5(i) = alpha_vector(3);
    alpha6(i) = alpha_vector(4);

    % solve for point acceleration
    aB(:, i) = FindAcc(a0, a, omega2, alpha2, e2, n2);
    aE(:, i) = FindAcc(a0, p, omega2, alpha2, eAE, nAE);
    aF(:, i) = FindAcc(a0, q, omega4(i), alpha4(i), eDF, nDF);
    aG(:, i) = FindAcc(aE(:, i), u, omega5(i), alpha5(i), e5, n5);

    % force analysis
    [eA2, nA2, LA2, s2A, s2B, s2E] = LinkCG(a, p, gamma2, xbar2, theta2(i));
    [eB3, nB3, LB3, s3B, s3C, ~] = LinkCG(b, 0, 0, xbar3, theta3(i));
    [eD4, nD4, LD4, s4D, s4C, s4F] = LinkCG(c, q, gamma4, xbar4, theta4(i));
    [eE5, nE5, LE5, s5E, s5G, ~] = LinkCG(u, 0, 0, xbar5, theta5(i));
    [eF6, nF6, LF6, s6F, s6G, s6W] = LinkCG(v, w, gamma6, xbar6, theta6(i));

    v2(:, i) = FindVel(v0, LA2, omega2, nA2);
    v3(:, i) = FindVel(vB(:, i), LB3, omega3(i), nB3);
    v4(:, i) = FindVel(v0, LD4, omega4(i), nD4);
    v5(:, i) = FindVel(vE(:, i), LE5, omega5(i), nE5);
    v6(:, i) = FindVel(vF(:, i), LF6, omega6(i), nF6);

    a2(:, i) = FindAcc(a0, LA2, omega2, alpha2, eA2, nA2);
    a3(:, i) = FindAcc(aB(:, i), LB3, omega3(i), alpha3(i), eB3, nB3);
    a4(:, i) = FindAcc(a0, LD4, omega4(i), alpha4(i), eD4, nD4);
    a5(:, i) = FindAcc(aE(:, i), LE5, omega5(i), alpha5(i), eE5, nE5);
    a6(:, i) = FindAcc(aF(:, i), LF6, omega6(i), alpha6(i), eF6, nF6);

    S_matrix = [U2    -U2    Z2    Z2    -U2    Z2    Z2    Z21;
                Z2     U2   -U2    Z2     Z2    Z2    Z2    Z21;
                Z2     Z2    U2   -U2     Z2    U2    Z2    Z21;
                Z2     Z2    Z2    Z2     U2    Z2   -U2    Z21;
                Z2     Z2    Z2    Z2     Z2   -U2    U2    Z21;
                s2A' -s2B'  Z12   Z12    -s2E' Z12   Z12      1;
                Z12   s3B' -s3C'  Z12    Z12   Z12   Z12      0;
                Z12   Z12   s4C' -s4D'   Z12   s4F'  Z12      0;
                Z12   Z12   Z12   Z12    s5E'  Z12  -s5G'     0;
                Z12   Z12   Z12   Z12    Z12  -s6F'  s6G'     0;];

    t_vector = [m2*a2(:, i);
                m3*a3(:, i);
                m4*a4(:, i);
                m5*a5(:, i);
                m6*a6(:, i) - FW;
                I2*alpha2;
                I3*alpha3(i);
                I4*alpha4(i);
                I5*alpha5(i);
                I6*alpha6(i) - dot(s6W, FW)];

    f_vector = S_matrix\t_vector;

    FA(:, i) = [f_vector(1); f_vector(2)];
    FB(:, i) = [f_vector(3); f_vector(4)];
    FC(:, i) = [f_vector(5); f_vector(6)];
    FD(:, i) = [f_vector(7); f_vector(8)];
    FE(:, i) = [f_vector(9); f_vector(10)];
    FF(:, i) = [f_vector(11); f_vector(12)];
    FG(:, i) = [f_vector(13); f_vector(14)];
    T2(i) = f_vector(15);

    % Energy method validation
    P_FW = dot(FW, vW(:, i));
    P_T2 = T2(i)*omega2;
    P_Ext(i) = P_FW + P_T2;

    P2 = InertialPower(m2, I2, v2(:, i), a2(:, i), omega2, alpha2);
    P3 = InertialPower(m3, I3, v3(:, i), a3(:, i), omega3(i), alpha3(i));
    P4 = InertialPower(m4, I4, v4(:, i), a4(:, i), omega4(i), alpha4(i));
    P5 = InertialPower(m5, I5, v5(:, i), a5(:, i), omega5(i), alpha5(i));
    P6 = InertialPower(m6, I6, v6(:, i), a6(:, i), omega6(i), alpha6(i));
    P_Kin(i) = P2 + P3 + P4 + P5 + P6;
end

% estimation of omega, velocity and accel values using numerical method
timestep = 2*pi/((N-1)*omega2);
omega3_estimated = FiniteDiffMethod(theta3, timestep);
omega4_estimated = FiniteDiffMethod(theta4, timestep);
vEx_estimated = FiniteDiffMethod(xE(1, :), timestep);
vEy_estimated = FiniteDiffMethod(xE(2, :), timestep);
vGx_estimated = FiniteDiffMethod(xG(1, :), timestep);
vGy_estimated = FiniteDiffMethod(xG(2, :), timestep);

aGx_estimated = FiniteDiffMethod(vG(1, :), timestep);
aGy_estimated = FiniteDiffMethod(vG(2, :), timestep);

% create a tiled plot
layout = tiledlayout(3, 2);

% plot velocity 
nexttile(2);
plot(theta2*180/pi, aG(1, :));
hold on
plot(theta2*180/pi, aGx_estimated, '.');
plot(theta2*180/pi, aG(2, :));
plot(theta2*180/pi, aGy_estimated, '.');
xlim([0 360]);
set(gca, 'xtick', 0:60:360);
xlabel('Crank Angle (degree)');
ylabel('Velocity (m/s^2)');
legend('analytical', 'estimated');

% plot validation from energy method
nexttile(3);
plot(theta2*180/pi, P_Kin);
hold on
plot(theta2*180/pi, P_Ext, '.');
xlim([0 360]);
legend('P\_Kin', 'P\_Ext', Location='best');
set(gca, 'xtick', 0:60:360);
xlabel('Crank Angle (degree)');
ylabel('Power (W)');

% plot torque
nexttile(4);
plot(theta2*180/pi, T2);
set(gca, 'xtick', 0:60:360);
xlim([0 360]);
xlabel('Crank Angle (degree)');
ylabel('Torque (Nm)');

%define linkage colour
cBlu = DefineColor([0, 110, 199]);
cBlk = DefineColor([0, 0, 0]);

% plot path
nexttile(1)
plot(xB(1, :), xB(2, :));
hold on
plot(xC(1, :), xC(2, :));
plot(xE(1, :), xE(2, :), 'Color', cBlu(1, :));
plot(xF(1, :), xF(2, :), 'Color', cBlu(7, :));
plot(xG(1, :), xG(2, :), 'Color', cBlk(5, :));
axis equal
xlim([-0.25, 0.25]);
ylim([-0.25, 0.25]);

%initialise linkage plot
if class == 2 || class == 1
    h1 = plot([x0(1, 1), xB(1, 1), xC(1, 1)], [x0(2, 1), xB(2, 1), xC(2, 1)]);
    h2 = plot([x0(1, 1), xE(1, 1), xB(1, 1)], [x0(2, 1), xE(2, 1), xB(2, 1)]);
    h3 = plot([xD(1, 1), xC(1, 1)], [xD(2, 1), xC(2, 1)]);
    h6 = plot([xD(1, 1), xF(1, 1), xC(1, 1)], [xD(2, 1), xF(2, 1), xC(2, 1)]);
    h7 = plot([xE(1, 1), xG(1, 1), xF(1, 1)], [xE(2, 1), xG(2, 1), xF(2, 1)]);
elseif class == 4
    h1 = plot([x0(1), xC(1, 1), xB(1, 1)], [x0(2), xC(2, 1), xB(2, 1)]);
    h2 = plot([xC(1, 1), xE(1, 1), xB(1, 1)], [xC(2, 1), xE(2, 1), xB(2, 1)]);
    h3 = plot([xD(1), xB(1, 1)], [xD(2), xB(2, 1)]);
end

h4 = patch([x0(1), xB(1, 1), xE(1, 1)], [x0(2), xB(2, 1), xE(2, 1)], cBlu(9, :), 'FaceAlpha', 0.5);
h8 = patch([xD(1), xC(1, 1), xF(1, 1)], [xD(2), xC(2, 1), xF(2, 1)], cBlu(10, :), 'FaceAlpha', 0.5);

%initialise joint marker
h5 = plot([x0(1), xB(1, 1), xC(1, 1), xD(1), xE(1, 1), xF(1, 1), xG(1, 1)], ...
    [x0(2), xB(2, 1), xC(2, 1), xD(2), xE(2, 1), xF(2, 1), xG(2, 1)], 'o', 'MarkerFaceColor', cBlk(5, :));

%initialise label
t1 = text(x0(1) + 0.005, x0(2), 0, "A");
t2 = text(xB(1, 1) + 0.005, xB(2, 1), 0, "B");
t3 = text(xC(1, 1) + 0.005, xC(2, 1), 0, "C");
t4 = text(xE(1, 1) + 0.005, xE(2, 1), 0, "E");
t5 = text(xD(1) + 0.005, xD(2), 0, "D");
t6 = text(xF(1, 1) + 0.005, xF(2, 1), 0, "F");
t7 = text(xG(1, 1) + 0.005, xG(2, 1), 0, "G");

% animation loop
j = 1;
direction = 1;
while j <= N
    if ~isgraphics(layout)
        break
    end

    if class == 2 || class == 1
        set(h1, 'XData', [x0(1), xB(1, j), xC(1, j)], 'YData', [x0(2), xB(2, j), xC(2, j)]);
        set(h2, 'XData', [x0(1), xE(1, j), xB(1, j)], 'YData', [x0(2), xE(2, j), xB(2, j)]);
        set(h3, 'XData', [xD(1), xC(1, j)], 'YData', [xD(2), xC(2, j)]);
        set(h6, 'XData', [xD(1), xF(1, j), xC(1, j)], 'YData', [xD(2), xF(2, j), xC(2, j)]);
        set(h7, 'XData', [xE(1, j), xG(1, j), xF(1, j)],'YData', [xE(2, j), xG(2, j), xF(2, j)]);
    elseif class == 4
        set(h1, 'XData', [x0(1), xC(1, j), xB(1, j)], 'YData', [x0(2), xC(2, j), xB(2, j)]);
        set(h2, 'XData', [xC(1, j), xE(1, j), xB(1, j)], 'YData', [xC(2, j), xE(2, j), xB(2, j)]);
        set(h3, 'XData', [xD(1), xB(1, j)], 'YData', [xD(2), xB(2, j)]);
    end

    %update joint marker
    set(h5, 'XData', [x0(1), xB(1, j), xC(1, j), xD(1), xE(1, j), xF(1, j), xG(1, j)], ...
    'YData', [x0(2), xB(2, j), xC(2, j), xD(2), xE(2, j), xF(2, j), xG(2, j)])

    %update coupler patch
    set(h4, 'XData', [x0(1), xB(1, j), xE(1, j)], 'YData', [x0(2), xB(2, j), xE(2, j)]);
    set(h8, 'XData', [xD(1), xC(1, j), xF(1, j)], 'YData', [xD(2), xC(2, j), xF(2, j)]);
    drawnow

    %update label position
    %check if text object still exist first, before update
    if all(isgraphics([t2, t3, t4]))
        set(t2, 'Position', [xB(1, j) + 0.005, xB(2, j)]);
        set(t3, 'Position', [xC(1, j) + 0.005, xC(2, j)]);
        set(t4, 'Position', [xE(1, j) + 0.005, xE(2, j)]);
        set(t6, 'Position', [xF(1, j) + 0.005, xF(2, j)]);
        set(t7, 'Position', [xG(1, j) + 0.005, xG(2, j)]);
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