function X_dot = getDot(X, t)
% This function uses a first order difference approximation to generate the
% time derivative state matrix X_dot
% INPUTS:
%   X - an (N x K) matrix representing the state
%   t - time vector
% OUTPITS:
%   X_dot - an (N x K-1) matrix representing the time derivative of the
%   state
    dt = t(2:end) - t(1:end-1);
    X_dot = (1./dt).*(X(:, 2:end) - X(:, 1:end-1));
end
