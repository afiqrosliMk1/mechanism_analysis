% Gear_fivebar_position_analysis.m
% conducts a position analysis on the geared fivebar linkage
% positions of points P and Q
% Afiq Rosli, 25/5/2025

% prepare workspace
clear variables; close all; clc;
a = 0.120;               % crank length (m)
b = 0.250;               % coupler 1 length (m)
c = 0.250;               % coupler 2 length (m)
d = 0.180;               % ground link length (m)
u = 0.150;               % length of link on gear 2 (m)
N1 = 24;                 % number of teeth on gear 1
N2 = 24;                 % number of teeth on gear 2
rho = N1/N2;             % gear ratio
phi = 0;                 % offset angle between gears
gamma3 = 20 * pi / 180;  % angle to point P on coupler 1
gamma4 = -20 * pi / 180; % angle to point Q on coupler 2
p = 0.200;               % distance to point P on coupler 1
q = 0.200;               % distance to point Q on coupler 2

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
xB = zeros(2, N);
xC = zeros(2, N);
xP = zeros(2, N);
xE = zeros(2, N);
xQ = zeros(2, N);

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

end

% plot path
fig = figure;
% plot(xB(1, :), xB(2, :));

plot(xQ(1, :), xQ(2, :));
hold on
plot(xP(1, :), xP(2, :));

%plot(xE(1, :), xE(2, :));

axis equal
xlim([-0.30, 0.40]);
ylim([-0.15, 0.25]);

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
    if ~ishandle(fig)
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

    %update label position
    set(t2, 'Position', [xB(1, j) + 0.005, xB(2, j)]);
    set(t3, 'Position', [xC(1, j) + 0.005, xC(2, j)]);
    set(t4, 'Position', [xE(1, j) + 0.005, xE(2, j)]);
    set(t6, 'Position', [xP(1, j) + 0.005, xP(2, j)]);
    set(t7, 'Position', [xQ(1, j) + 0.005, xQ(2, j)]);

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