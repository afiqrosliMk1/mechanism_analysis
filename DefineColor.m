% DefineColor.m
% makes a palette of custom colors for plotting in MATLAB
% the rist color in the palette gives the full base color, and the rest
% fades gradually to white
%
% colorBase = 1x3 input array giving RGB values (between 0 and 255)
% C = 10x3 output array giving RGB values (between 0 and 1)

function C = DefineColor(colorBase)
    C = zeros(11, 3);
    for i = 1:11
        for j = 1:3
            C(i, j) = (255 - colorBase(j)) / 10 * (i - 1) + colorBase(j);
        end
    end

    C = C/255;
end