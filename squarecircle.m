% squarecircle.m
% makes a plot of eight squares arranged in a circle

clear variables; close all; clc;

unitSquare = 0.5 * [ -1 1 1 -1 -1;
                     -1 -1 1 1 -1];

for i = 1:8
    % angle at which to place square, i - 1 to ensure it starts at 0
    theta = (i - 1) * 2 * pi / 8; 
    
    % 2 is the length of vector r. the two values are components of vector
    r = 2 * [cos(theta); sin(theta)]; 

    % create new square center at end of vector r
    newSquare = unitSquare + r;

    %fill and edge color of square. [1 1 1] - white. [0 0 0] - black
    % rgb color vector is multiplied by (i -1)/8
    col = (i - 1) * [ 1 1 1] / 8;

    fill(newSquare(1, :), newSquare(2, :), col, 'EdgeColor', col, 'LineWidth', 2);
    hold on
end

xlabel('x (m)'); xlim([-3 3])
ylabel('y (m)'); ylim([-3 3])
grid on
axis equal