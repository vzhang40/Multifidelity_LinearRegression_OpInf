function x0 = getIC(x, param)
% This function 
% INPUTS: 
%   x - (N x 1) discretized spatial coordinate
%   param - (2 x pq) parameters [a; b] assoicated with each intial
%       condition
% OUTPUTS:
%   x0 - (N x 1) discretized initial condition based on parameters
    x01 = param(1,:)'*sin(2.*pi.*x./x(end));
    x02 = cos(2*pi*param(2, :)'*x./x(end));
    x0 = (x01 + x02)';
end
