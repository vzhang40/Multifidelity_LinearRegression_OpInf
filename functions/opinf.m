function Ohat = opinf(Xr, Xrdot, m, XXT)
% This function returns the learned reduced operator(s)
%
% INPUTS: 
%   XXT: (r x r) inverse term in estimating operators/parameters
%   XXdot: (r x r) derivative data, can be (1 x K) if estimating operators
%   by row
%   K: number of samples used
% where r is the reduced dimension of concatenated operators

% OUTPUTS:
%   Ohat: (r x r) reduced operators(s) can be (1 x r) if estimating
%   operators by row
p = size(XXT, 1);
r = size(Xrdot, 1);
Ohat = zeros(r, p);
for i = 1:r
    % Single fi Data
    indsSF = randperm(size(Xr, 2), m);
    X = Xr(:, indsSF);
    Y = Xrdot(:, indsSF)';
    beta_single = (XXT+1e-3*eye(p))\(1./size(Y, 1).*X*Y(:, i).*1);
    Ohat(i, :) = beta_single';
end
end