% Fourbar_Position_Analysis_with_non_grashoff.m
% ++with solver for non grashoff case
% conducts a position analysis of the fourbar linkage and plots the
% positions of points B, C and P
% Afiq Rosli, 10/5/2025
% the methods of finding coupler angle (theta 3) and rocker angle (theta 4)
% is taken from Karl B Dyer et al

% prepare workspace
clear variables; close all; clc;
a = 0.17;               % crank length (m)
b = 0.2;               % coupler length (m)
c = 0.13;               % rocker length (m)
d = 0.22;              % ground link length (m)
p = 0.150;               % length from B to P (m)
gamma = 20 * pi / 180; % angle between BP and coupler (radian)

% if grashoff condition met, grashoff = 0
% if non-grashoff, grashoff = 1
grashoff = 0;
%ccw = 0, cw = 1;
cw = 0;
% check for grashoff condition
[S, L, P, Q] = slpq(a , b, c, d);
if S + L < P + Q 
    fprintf("grashoff\n")
    theta2_min = 0;
    theta2_max = 2*pi;
else
    if L == a
        fprintf("non grashoff class 7 - longest link is crank (link ""a"")\n")
        grashoff = 2;
    elseif L == b
        fprintf("non grashoff class 8&9 - longest link is coupler (link ""b"")\n")
        grashoff = 3;
        cw = 1;
        theta2_max = acos((a^2 + d^2 - (b - c)^2) / (2*a*d));
        theta2_min = -theta2_max;
    elseif L == c
        fprintf("non grashoff class 10 - longest link is rocker (link ""c"")\n")
        grashoff = 4;
    elseif L == d
        fprintf("non grashoff class 5&6 - longest link is ground (link ""d"")\n")
        grashoff = 1;
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
xD = [-d; 0];
xB = zeros(2, N);
xC = zeros(2, N);
xP = zeros(2, N);

% main loop
for i = 1:N
    if cw == 0
        theta2(i) = (i - 1) * (theta2_max - theta2_min) / (N - 1) + theta2_min;
    elseif cw == 1
        theta2(i) = -1 * (i - 1) * (2*pi - (theta2_max - theta2_min)) / (N - 1) + theta2_min;
    end
  
    %r = d - a * cos(theta2(i));
    r = d + c * cos(theta2(i));
    %s = a * sin(theta2(i));
    s = c * sin(theta2(i));
    
    f = sqrt(r^2 + s^2);
    
    % open configuration
    %cos_val = (b^2 + c^2 - f^2) / (2*b*c);
    cos_val = (a^2 + b^2 - f^2) / (2*a*b);
    % Clamp to [-1, 1]. Important to handle floating point error. otherwise
    % delta will give you 0 + 2.107342425544702e-08i for theta2_min angle.
    cos_val = min(1, max(-1, cos_val));   
    delta = acos(cos_val);
    
    % g = b - c * cos(delta);
    % h = c * sin(delta);
    g = b - a * cos(delta);
    h = a * sin(delta);
    
    %changing the "x" to negative somehow works. but why?
    theta3_prime(i) = atan2((h*r - g*s), (g*r + h*s));
    theta3(i) = pi - theta3_prime(i);
    
    %theta4(i) = theta3(i) + delta;
    theta4(i) = theta3(i) - delta;
    
    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));
    [e4, n4] = UnitVector(theta4(i));
    %[eBP, nBP] = UnitVector(theta3(i) + gamma);
    [eCP, nCP] = UnitVector(theta3(i) - gamma);
    
    %xB(:, i) = FindPos(x0, c, e2);
    xB(:, i) = FindPos(xD, a, e4);
    %xC(:, i) = FindPos(xD, a, e4);
    xC(:, i) = FindPos(x0, c, e2);
    %xP(:, i) = FindPos(xB(:, i), p, eBP);
    xP(:, i) = FindPos(xC(:, i), p, eCP);
    

end

% plot path
fig = figure;
plot(xB(1, :), xB(2, :));
hold on
plot(xC(1, :), xC(2, :));
plot(xP(1, :), xP(2, :));
axis equal

%initialise linkage plot
%h1 = plot([x0(1), xB(1, 1), xC(1, 1)], [x0(2), xB(2, 1), xC(2, 1)]);
h1 = plot([x0(1), xC(1, 1), xB(1, 1)], [x0(2), xC(2, 1), xB(2, 1)]);
%h2 = plot([xB(1, 1), xP(1, 1), xC(1, 1)], [xB(2, 1), xP(2, 1), xC(2, 1)]);
h2 = plot([xC(1, 1), xP(1, 1), xB(1, 1)], [xC(2, 1), xP(2, 1), xB(2, 1)]);
%h3 = plot([xD(1), xC(1, 1)], [xD(2), xC(2, 1)]);
h3 = plot([xD(1), xB(1, 1)], [xD(2), xB(2, 1)]);
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
    %set(h1, 'XData', [x0(1), xB(1, j), xC(1, j)], 'YData', [x0(2), xB(2, j), xC(2, j)]);
    set(h1, 'XData', [x0(1), xC(1, j), xB(1, j)], 'YData', [x0(2), xC(2, j), xB(2, j)]);
    %set(h2, 'XData', [xB(1, j), xP(1, j), xC(1, j)], 'YData', [xB(2, j), xP(2, j), xC(2, j)]);
    set(h2, 'XData', [xC(1, j), xP(1, j), xB(1, j)], 'YData', [xC(2, j), xP(2, j), xB(2, j)]);
    %set(h3, 'XData', [xD(1), xC(1, j)], 'YData', [xD(2), xC(2, j)]);
    set(h3, 'XData', [xD(1), xB(1, j)], 'YData', [xD(2), xB(2, j)]);


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
    if grashoff == 1
        if (j >= N && direction > 0) || (j <= 1 && direction < 0)
            direction = direction * -1;
        end
    elseif grashoff == 2
        fprintf("under construction\n")
    elseif grashoff == 3
        if (j >= N && direction > 0) || (j <= 1 && direction < 0)
            direction = direction * -1;
        end
    elseif grashoff == 4
        fprintf("under construction\n")
    elseif grashoff == 0
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