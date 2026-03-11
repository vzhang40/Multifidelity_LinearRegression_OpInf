function [m1, m2] = mfmcOptimalSampleSize(p, w, Xhrdot, Xlrdot)
% This function takes the state and state derivative matrices as well as
% budget and model costs to determine the mfmc optimal sample allocation
%
% INPUTS: 
%   p - (1 x P) vector of budgets
%   w - vector of model costs
%   Xhrdpt - (r x K) high fidelity state derivative matrix
%   Xlrdot - (r x K) low fidelity state derivative matrix
%       where K is the number of samples to estimate statistics
% 
% OUTPUTS:
%   m1: (r x P) matrix of optimal high fidelity sample sizes 
%   m2: (r x P) matrix of optimal low fidelity sample sizes 
r = size(Xhrdot, 1);
rhos = zeros(r, 3);
m1 = zeros(r, length(p));
m2 = zeros(r, length(p));
rhos(:, 1) = ones(r, 1);
for i = 1:r
    temp = corrcoef(Xhrdot(i, :), Xlrdot(i, :));
    rhos(i, 1:end-1) = temp(1, :);
    rs = sqrt( w(1)*(rhos(i, 1:end-1).^2 - rhos(i, 2:end).^2) ./ (w'*(1 - rhos(i, 2).^2)))';
    m1(i, :) = floor(p./(w'*rs));
    m2(i, :) = floor(m1(i, :).*rs(2));
end
end