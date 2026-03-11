function alpha = optimalControlVariate(bigG)
% This function takes the state and state derivative matrices in order to
% find optimal scalar control variates term to estimate each row of a
% reduced order operator of input size r.
%
% INPUTS: 
%   bigG - (2 x 2 x r) cell of covariance matrices
% OUTPUTS:
%   alpha - (r x 1) vector of optimal scalar control variates for latent
%   dimension of derivative term

r = size(bigG, 3);
alpha = zeros(r, 1);
for i = 1:r
    alpha(i) = trace(bigG{1, 2, i}) / trace(bigG{2, 2, i});
end
end
