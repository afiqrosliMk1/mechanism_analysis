% squareplot.m
% plots a filled square

clear variables; close all; clc

unitSquare = 0.5 * [-1 1 1 -1 -1;
                    -1 -1 1 1 -1];

%plot(unitSquare(1, :), unitSquare(2, :), 'LineWidth', 2)
fill(unitSquare(1, :), unitSquare(2, :), 'b', 'LineWidth', 2);
hold on

xlabel('x (m)'); xlim([-2, 2]);
ylabel('y (m)'); ylim([-2, 2]);

grid on
axis equal

newSquare = 2 * unitSquare + [1; 1];
fill(newSquare(1, :), newSquare(2, :), 'r', 'LineWidth', 2, 'FaceAlpha', 0.5);