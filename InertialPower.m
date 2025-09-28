% InertialPower.m
% computes the total intertial power of a link
% m = mass of link
% I = moment of intertia of link
% v = velocity of center of mass of link
% a = acceleration of center of mass of link
% omega = angular velocity of link
% alpha = angular acceleration of link

function P = InertialPower(m, I, v, a, omega, alpha)
    P = m*dot(v, a) + I*omega*alpha; %intertial power of link
end

%how to use it?
%for each link, compute the external power supplied by external forces and
%moment
% PF = dot(FP, vP); %power from external force FP 
% PT = T2 * omega2; %power from driving torque T2
% PExt = PF + PT; %total external power
% 
% %then calculate the inertial power
% P2 = InertialPower(m2, I2, v2, a2, omega2, alpha2);
% P3 = InertialPower(m3, I3, v3, a3, omega3, alpha3);
% P_kin = P2 + P3; %total intertial power
% 
% %then make a plot to overlay both plots. Energy given to links, PExt must be
% %exactly as energy absorbed by the links, P_kin.
% plot(theta*180/pi, P_kin);
% hold on
% plot(theta*180/pi, PExt);