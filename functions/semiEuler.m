function X = semiEuler(A, F, x0, time)
% This function uses the semi-implicit euler's integration scheme.
% INPUTS:
%   A - (N x N) linear operator
%   F - (N x (N(N+1)/2) quadratic operator
%   x0 - (N x 1) initial condition
%   time - (1 x K) vector of times
% OUTPUTS:
%   X - (N x K) integrated state matrix

N = size(A, 1);
K = length(time);
X = zeros(N, K);
X(:, 1) = x0; 

% iterate through each time step
for i = 1:K-1
    dt = time(i+1) - time(i);
    new_X = (eye(N) - dt.*A)\(X(:, i) + dt*F*uniKron(X(:, i)));
    X(:, i+1) = new_X;
end
end
