% threebar_position_analysis.m
% Conducts a position analysis on the threebar crank-slider linkage
% by Afiq Rosli, 17/11/2024

% Prepare Workspace
clear variables; close all; clc;

% Linkage dimensions
a = 0.100; % crank length (m)
d = 0.150; % length between ground pins (m)
p = 0.300; % slider length (m)

% Ground pins
x0 = [0; 0]; % point A (the origin) - vector
xD = [d; 0]; % point D - vector

N = 361; % number of times to perform position calculations

theta2 = zeros(1, N); % allocate space for crank angle
theta3 = zeros(1, N); % allocate space for crank angle
b = zeros(1, N); % allocate space for slider length

xB = zeros(2, N);
xP = zeros(2, N);

% Main Loop
for i = 1: N
    theta2(i) = ( i - 1 ) * (2 * pi) / ( N - 1); % (2*pi)/(N-1) gives you angle per step. 
end
