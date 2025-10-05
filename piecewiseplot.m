clc; close all; clear variables;
% demo how to plot piecewise
% especially useful for acceleration profile that has initial pulse

T = 2; 
lambda = 2*pi/T;
N = 1000;
dt = T/(N-1);

t = 0: dt: 3;
theta2 = zeros(size(t));

theta2(t>=0 & t<=T) = 0.25*(lambda*t(t>=0 & t<=T) - sin(lambda*t(t>=0 & t<=T)));
theta2(t>T) = pi/2;
omega2(t>=0 & t<=T) = 0.25*lambda*(1 - cos(lambda*t(t>=0 & t<=T)));
omega2(t>T) = 0;
alpha2(t>=0 & t<=T) = 0.25*lambda^2*sin(lambda*t(t>=0 & t<=T));
alpha2(t>T) = 0;

layout = tiledlayout(3, 1);

nexttile;
plot(t, theta2*180/pi);
xlim("auto");
ylim("auto");
set(gca, 'ytick', 0:45:90)
xlabel('time (s)');
ylabel('crank angle (degree)');

nexttile;
plot(t, omega2);
xlim("auto");
ylim("auto");
xlabel('time (s)');
ylabel('angular velocity (rad/s)');

nexttile;
plot(t, alpha2);
xlim("auto");
ylim("auto");
xlabel('time (s)');
ylabel('angular acceleration (rad/s^2)');