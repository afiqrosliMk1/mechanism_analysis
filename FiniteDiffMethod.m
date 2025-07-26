% FINITEDIFFMETHOD gives derivative estimate using forward difference method
% f = function/array of data points whose derivatives is to be estimated
% derivative = the calculated derivative
% dt = time step
function [derivative] = FiniteDiffMethod(f, delta_t)
    N = length(f);
    derivative = zeros(1, N);

    for i = 1: N - 1
        derivative(i) = (f(i + 1) - f(i))/delta_t;
    end

    % assign the final value of derivative equal to previous one, because
    % we had reached the end
    derivative(N) = derivative(N - 1);
end