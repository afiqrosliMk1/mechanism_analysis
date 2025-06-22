% sixbar_Position_Analysis.m
% conducts a position analysis of the sixbar linkage and plots the
% positions of points E, F and G
% Afiq Rosli, 21/6/2025
% the methods of finding coupler angle (theta 3) and rocker angle (theta 4)
% is taken from Karl B Dyer et al

% prepare workspace
clear variables; close all; clc;
% crank length (m)
a = 0.070;
% coupler length (m)
b = 0.100;
% rocker length (m)
c = 0.09;
% ground link length (m)
d = 0.110;
% length from B to E (m)
p = 0.150; 
% length from D to F (m)
q = 0.150;
% length from E to G (m)
u = 0.120;
% length from G to F (m)
v = 0.160;
% angle between BE and a (radian)
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
xE = zeros(2, N);
xF = zeros(2, N);
xG = zeros(2, N);
alpha = 0;
beta = 0;
theta2_prime = 0;
theta5_prime = 0;
theta6_prime = 0;
theta5 = zeros(1, N);
theta6 = zeros(1, N);

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
        [eAE, nAE] = UnitVector(theta2(i) + gamma);
        [eDF, nDF] = UnitVector(theta4(i) - gamma);
    elseif class == 4
        [eCP, nCP] = UnitVector(theta2(i) - gamma);
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
end

%define linkage colour
cBlu = DefineColor([0, 110, 199]);
cBlk = DefineColor([0, 0, 0]);

% plot path
fig = figure;
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
    if ~isgraphics(fig)
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