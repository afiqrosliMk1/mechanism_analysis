% LinkCG.m
% calculates the vectors associated with the center of mass of two or three
% pin link
%
% **inputs**
% a = length of link from pin A to pin B
% c = length of link from pin A to pin C (zero if two-pin link)
% gamma = angle between AB and AC (zero if two-pin link)
% xbar = 2x1 vector that contains [xbar, ybar] - coordinate of center of
% mass in LOCAL coordinate system
% theta = angle of link relative to global coordinate system. measured from
% global x-axis to AB
%
% **outputs**
% eAg, nAg = unit vector and unit normal FROM point A TO center of mass
% LAg = length from point A to center of mass
% sgA = normal vector to vector rgA
% sgB = normal vector to vector rgB
% sgC = normal vector to vector rgC
%
% **others to note**
% rgA = vector FROM center of mass TO point A
% rgB = vector FROM center of mass TO point B
% rgC = vector FROM center of mass TO point C
% suffix _prime = meaning that that particular vector is relative to LOCAL
% Coordinate sytem which x-axis is alinged with AB and y-axis point up
%
% Note: vector in local coordinate is eventually tranformed into global
% coordinate by multiplying with Rotation Matrix

function [eAg, nAg, LAg, sgA, sgB, sgC] = LinkCG(a, c, gamma, xbar, theta)
    rAB_prime = a * [1; 0];
    rAC_prime = c * [cos(gamma); sin(gamma)];
    rgA_prime = -[xbar(1); xbar(2)];
    rgB_prime = rgA_prime + rAB_prime;
    rgC_prime = rgA_prime + rAC_prime;

    sgA_prime = [-rgA_prime(2); rgA_prime(1)];
    sgB_prime = [-rgB_prime(2); rgB_prime(1)];
    sgC_prime = [-rgC_prime(2); rgC_prime(1)];

    LAg = norm(rgA_prime);
    % handle situation where point coincide with center of mass, LAg = 0,
    % so eAg points to nowhere [0; 0]. we use 1e-12 (very small number) 
    % instead of zero to handle floating points rounding error.
    if LAg < 1e-12
        eAg_prime = [0; 0];
        nAg_prime = [0; 0];
    else
        eAg_prime = -(1/LAg)*rgA_prime;
        nAg_prime = [-eAg_prime(2); eAg_prime(1)];
    end

    % Rotation matrix
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

    % Transform to global coordinate system
    sgA = R * sgA_prime;
    sgB = R * sgB_prime;
    sgC = R * sgC_prime;
    eAg = R * eAg_prime;
    nAg = R * nAg_prime;
end