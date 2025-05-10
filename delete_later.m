clear variables; close all; clc;

N = 3;
j = 1;
direction = 1;
while j <= N
    fprintf("%d\n", j);
    if (j >= N && direction > 0) || (j <= 1 && direction < 0)
        direction = direction * -1;
    end
    j = j + 1*direction;
end