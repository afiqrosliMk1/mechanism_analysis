clear variables; close all; clc
red_shade = [1, 0, 0]

for i = 1:1:8
    %remember from best practice, x and y coordinate we manage it in a
    %matrix
    square = [1, -1, -1, 1, 1;
        1, 1, -1, -1, 1];
    radius = 5;
    %i got the pattern below by writing on a piece of paper and see how i
    %is related to each increment of angle theta
    theta = pi * (i - 1) / 4;

    %fprintf("%.6f\n", angle)
    x_offset = radius * cos(theta);
    y_offset = radius * sin(theta);
    %r is the vector from origin to the center of square
    r = [x_offset; 
        y_offset];
    square = square + r;
    plot(square(1,:), square(2,:))

    %change the shade of red according to each loop index, i.
    %max rgb value is 1. eg: [1, 1, 1]
    %red_shade(1) = red_shade(1) - red_shade(1) / 8
    fprintf("%d", i)
    red_shade = red_shade - (i / 32) * red_shade
    fill(square(1,:), square(2,:), red_shade)
    hold on
end

xlim([-10, 10])
ylim([-10, 10])
%make axis aspect ratio equal
axis equal
hold off