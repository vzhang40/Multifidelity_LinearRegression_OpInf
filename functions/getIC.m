function x0 = getIC(x, xi)
% This function 
% INPUTS: 

% OUTPUTS:
    s = 3/2;
    n = 1:(size(xi, 1)/2);
    a = (x(end)./(2*n*pi)).^s;
    X = [sin(2*pi*n.*x'./x(end)).*a, cos(2*pi*n.*x'./x(end)).*a];    % this is a sin basis multiplied by decaying weights
    x0 = sqrt(2)*X*xi; 
end
