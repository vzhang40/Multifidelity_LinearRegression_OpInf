function operators = getMatrices(N, mu, flag)
% This function gets the operator matrices for PDE examples: 
% 1) 1D Periodic Heat Equation (Linear)
% 2) 1D Periodic Viscous Burger's Equation (Quadratic)
%
% INPUTS: 
%   N: spatial dimension size of desired operator (length of state x)
%   mu: pde parameter (heat diffusivity, viscocity) 
%   flag: pde identifier (heat or Burger)
% OUTPUTS: 
%   operators: structure with possible fields:
%       A: (N x N) Linear Operator
%       B: (N x w) Control Operator
%       F: (N x (n+1)C(2)) Quadratic Operator in Unique Kronecker Product

dx = 1/(N + 1);
A = (mu/dx^2)*(diag(-2*ones(N, 1), 0) + diag(ones(N-1, 1), -1) + diag(ones(N-1, 1), 1) + diag(1, -(N-1)) + diag(1, N-1));
operators.A = A;

if flag == "Burger"
    H = zeros(N, N^2);

    for i = 1:N
    index = N*(i-1) + i;
    if i ~= 1
        H(i, index - 1) = -1;
    else
        H(i, end - (N-1)) = -1;
    end
    
    if i ~= N
        H(i, index + 1) = 1;
    else
        H(i, N) = 1;
    end
    end
    
    H = (1/(2*dx)).*H;
    operators.H = H;
    operators.F = eliminate(H);
end
 
end
