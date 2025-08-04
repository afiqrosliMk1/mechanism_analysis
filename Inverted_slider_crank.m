% Inverted_slider_crank.m
% conducts a position analysis of the slider crank and plots the
% positions of points B, C and P
% Afiq Rosli, 25/5/2025

% prepare workspace
clear variables; close all; clc;
a = 0.080;               % crank length (m)
c = 0.130;               % rocker length (m)
d = 0.200;              % ground link length (m)
p = 0.350;               % length from B to P (m)
delta = 60 * pi/ 180;  % angle between slider and rocker (radian)
h = c * sin(delta); % h is constant. only calculate it once
eps = 0.000001; % tiny number added to theta2 to keep its bounds

% if grashoff condition met, grashoff = true
if ((a + d) < c * sin(delta))
    fprintf("linkage cannot be assembled\n")
else
    if abs(a - d) >= c * sin(delta)
        grashoff = true;
        fprintf("linkage is grashoff\n")
        theta2_min = 0;
        theta2_max = 2*pi;
    else
        grashoff = false;
        fprintf("non grashoff\n")
        theta2_min = acos((a^2 + d^2 - (c*sin(delta))^2) / (2*a*d))  + eps;
        theta2_max = 2*pi - theta2_min;
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
x0 = [0; 0];
xD = [d; 0];
xB = zeros(2, N);
xC = zeros(2, N);
xP = zeros(2, N);

%variable for velocity analysis
omega2 = 10; %rad/s
omega3 = zeros(1, N);
bdot = zeros(1, N);
v0 = [0; 0]; %velocity of pin A
vB = zeros(2, N);
vP = zeros(2, N);

% main loop
for i = 1:N
    theta2(i) = (i - 1) * (theta2_max - theta2_min) / (N - 1) + theta2_min;
  
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
    
    xB(:, i) = FindPos(x0, a, e2);
    xP(:, i) = FindPos(xB(:, i), p, e3);
    xC(:, i) = FindPos(xD, c, e4);

    % velocity analysis
    A_matrix = [(b(i)*n3 - c*n4), e3];
    b_vector = -a*omega2*n2;
    omega_vector = A_matrix\b_vector;
    omega3(i) = omega_vector(1);
    bdot(i) = omega_vector(2);

    %calculate point velocity
    vB(:, i) = FindVel(v0, a, omega2, n2);
    vP(:, i) = FindVel(vB(:, i), p, omega3(i), n3);

end

%--lines of code below are purely for validation of above point
%velocity. using finite difference method
timestep = 2*pi/((N-1)*omega2);
vPx_estimated = FiniteDiffMethod(xP(1, :), timestep);
vPy_estimated = FiniteDiffMethod(xP(2, :), timestep);

% plot velocity
fig2 = figure;
plot(theta2*180/pi, vP(1, :));
hold on
%plot estimated value from finite difference
plot(theta2*180/pi, vPx_estimated(1, :), '.');
xlim([0 360]);
set(gca, 'xtick', 0:60:360)
set(fig2, 'position', [50, 200, 600, 400])
xlabel('crank angle (degree)');
legend('vPx', 'vPx\_estimated')

% plot path
fig1 = figure;
plot(xB(1, :), xB(2, :));
hold on
plot(xC(1, :), xC(2, :));
plot(xP(1, :), xP(2, :));
axis equal
set(fig1, 'position', [700, 200, 600, 400])

%initialise linkage plot
h1 = plot([x0(1), xB(1, 1), xC(1, 1), xP(1, 1)], [x0(2), xB(2, 1), xC(2, 1), xP(2, 1)]);
h3 = plot([xD(1), xC(1, 1)], [xD(2), xC(2, 1)]);

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
    if ~isgraphics(fig1)
        break
    end
    set(h1, 'XData', [x0(1), xB(1, j), xC(1, j), xP(1, j)], 'YData', [x0(2), xB(2, j), xC(2, j), xP(2, j)]);
    set(h3, 'XData', [xD(1), xC(1, j)], 'YData', [xD(2), xC(2, j)]);

    %update joint marker
    set(h5, 'XData', [x0(1), xB(1, j), xC(1, j), xD(1), xP(1, j)], ...
    'YData', [x0(2), xB(2, j), xC(2, j), xD(2), xP(2, j)])

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