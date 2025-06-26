% Fourbar_using_newton_raphson.m
% solves the fourbar linkage using Newton-Raphson algorithm
% Afiq Rosli 26/6/2025

% prepare workspace
clear variables; close all; clc

% linkage dimension
% crank length (m)
a = 0.130;
% coupler length (m)
b = 0.200;
% rocker length (m)
c = 0.170;
% ground (m)
d = 0.220;
% length from B to P (m)
p = 0.150;
% angle between BP and coupler (rad)
gamma = 20*pi/180;

% ground pins
% ground pin at A
x0 = [0; 0]; 
% ground pin at D
xD = [d; 0];

% number of times to perform calculation
N = 361;
% point B, C and P
xB = zeros(2, N);
xC = zeros(2, N);
xP = zeros(2, N);

theta2 = zeros(1, N);
theta3 = zeros(1, N);
theta4 = zeros(1, N);

% initial guess
t3 = pi/4;
t4 = pi/2;

for i = 1:N
    theta2(i) = (i - 1)*(2*pi)/(N - 1);

    % Newton-Raphson method
    for j = 1:5
        phi(1, 1) = a*cos(theta2(i)) + b*cos(t3)  - c*cos(t4) - d;
        phi(2, 1) = a*sin(theta2(i)) + b*sin(t3) - c*sin(t4);

        % if constraint equations are close enough to zero
        if (norm(phi) < 0.00001)
            theta3(i) = t3;
            theta4(i) = t4;
            break;
        end

        % calculate Jacobian matrix
        J = [-b*sin(t3), c*sin(t4); 
              b*cos(t3), -c*cos(t4)];

        % update variables using Newton-Raphson equation
        delta_t = -J\phi;
        t3 = t3 + delta_t(1);
        t4 = t4 + delta_t(2);
    end

    % calculate unit vectors
    [e2, n2] = UnitVector(theta2(i));
    [e3, n3] = UnitVector(theta3(i));
    [e4, n4] = UnitVector(theta4(i));
    [eBP, nBP] = UnitVector(theta3(i) + gamma);

    % solve for position B, C and P
    xB(:, i) = FindPos(x0, a, e2);
    xC(:, i) = FindPos(xD, c, e4);
    xP(:, i) = FindPos(xB(:, i), p, eBP);
end

plot(xB(1, :), xB(2, :));
hold on
plot(xC(1, :), xC(2, :));
plot(xP(1, :), xP(2, :));
axis equal
title('paths of point B, C and P on 4 bar linkage')
legend('Point B', 'Point C', 'Point P', 'Location', 'southeast')