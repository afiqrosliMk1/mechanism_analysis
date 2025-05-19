% Fourbar_Position_Analysis_with_non_grashoff.m
% ++with solver for non grashoff case
% conducts a position analysis of the fourbar linkage and plots the
% positions of points B, C and P
% Afiq Rosli, 10/5/2025
% the methods of finding coupler angle (theta 3) and rocker angle (theta 4)
% is taken from Karl B Dyer et al

% prepare workspace
clear variables; close all; clc;
% crank length (m)
%a = 0.13;
a = 0.22;
% coupler length (m)
b = 0.2;
% rocker length (m)
c = 0.17;
% ground link length (m)
%d = 0.22;
d = 0.13;
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

end

% plot path
fig = figure;
plot(xB(1, :), xB(2, :));
hold on
plot(xC(1, :), xC(2, :));
plot(xP(1, :), xP(2, :));
axis equal

%initialise linkage plot
if class == 2 || class == 1
    h1 = plot([x0(1), xB(1, 1), xC(1, 1)], [x0(2), xB(2, 1), xC(2, 1)]);
    h2 = plot([xB(1, 1), xP(1, 1), xC(1, 1)], [xB(2, 1), xP(2, 1), xC(2, 1)]);
    h3 = plot([xD(1), xC(1, 1)], [xD(2), xC(2, 1)]);
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
    if ~isgraphics(fig)
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