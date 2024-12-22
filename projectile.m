% projectile.m
% animates the path of a projectile

clear variables; close all; clc

t = 0:0.1:15; % [s] vector of times
vx0 = 50;     % [m/s] initial velocity in the x-direction
vy0 = 50;     % [m/s] initial velocity in the y-direction
g = 9.81;     % [m/s^2] acceleration of gravity

% open plot window without placing anything yet
figure

set(gcf, 'Position', [50 50 1200 500])

for i = 1:length(t)
    x = vx0 * t(i);
    y = vy0 * t(i) - 0.5 * g * t(i) ^ 2;
    plot(x, y, 'o', 'MarkerFaceColor','b')

    axis equal
    xlabel('x (m)'); xlim([0 600])
    ylabel('y (m)'); ylim([0 150])
    grid on

    % enable drawnow to see animation, enable hold on/off to see static
    % plot
    drawnow
    hold on
end
%hold off