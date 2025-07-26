% Function FindVel.m
% calculates the translational velocity at a point on the linkage
% using relative velocity formula
%
% v0 = velocity of the first point
% L = length of vector to second point on the link
% omega = angular velocity of the link
% n = unit normal to vector between first and second point
% v = velocity of the second point

function v = FindVel(v0, L, omega, n)
    v = v0 + L* omega * n;
end