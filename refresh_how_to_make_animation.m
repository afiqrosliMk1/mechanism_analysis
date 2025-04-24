clear variables; close all; clc;

t = 0:0.1:10;
v0x = 50;
v0y = 50;
g = 9.81;

figure
set(gcf, 'Position', [50 50 1200 500])

for i = 1:length(t)
    x = v0x * t(i);
    y = v0y * t(i) - 0.5 * g * t(i)^2;

    plot(x, y, 'o', MarkerFaceColor='g');
    axis equal
    %you need to define xlim and ylim, otherwise you won't get a static
    %frame for your animation, and it won't work 
    xlim([0 600])
    ylim([0 150])

    drawnow
    
end
