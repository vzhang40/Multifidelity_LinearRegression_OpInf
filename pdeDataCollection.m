%% This script collects 120000 snapshot vectors to 
%
%

clear; close all; clc;

%% Model Dimensions and Parameters
L = 1; % length
T = 1.5; % final time
deta = [2^-7, 2^-4]; % mesh width for discretization
dt = [10^-3, 10^-1]; % time step size
mu = 0.01; % thermal diffusivity
N = L./deta; % number of spatial grid points
K = T./dt; % number of time grid points

%% Dimension Set Up
timesteps = 1; 

x = linspace(0, L, N(1));
xl = linspace(0, L, N(2));
t = linspace(0, T, K(1)+1);
tl = linspace(0, T, K(2)+1);

t = t(1:timesteps + 1);
tl = tl(1:timesteps + 1);

p = 6000; q = 20;
para = getRandIC(p, q);
x0 = getIC(x, para);
x0l = getIC(xl, para);
x0l = interp1(xl, x0l, x, "linear"); 

datah = zeros(N(1), p*q);
datahdot = zeros(N(1), p*q);
datal = zeros(N(1), p*q);
dataldot = zeros(N(1), p*q);

%% Generating HEAT EQ Data
operators = getMatrices(N(1), mu, "heat");

disp("Start Heat Data Collection")
for i = 1:p*q

% High fidelity data
Xhtemp = backwardEuler(operators.A, x0(:, i), t);
Xhdottemp = getDot(Xhtemp, t);
Xhtemp = Xhtemp(:, 2:end);
datah(:, i) = Xhtemp;
datahdot(:, i) = Xhdottemp;

% Low fidelity data
Xltemp = backwardEuler(operators.A, x0l(:, i), tl);
Xltemp = interp1(tl, Xltemp', t)';
Xldottemp = getDot(Xltemp, t);
Xltemp = Xltemp(:, 2:end);
datal(:, i) = Xltemp;
dataldot(:, i) = Xldottemp;

% Tracking
if mod(i, p*q./20) == 0
    disp("Heat Percent Done: " + num2str(i*100./(p*q)) + "%")
end

end

save("heatEqData.mat","operators","datah", "datahdot", "datal", "dataldot", "para")
disp("Variables saved to heatEqData.mat")

%% Generating BURGER EQ Data

% initializing
datah = zeros(N(1), p*q);
datahdot = zeros(N(1), p*q);
datal = zeros(N(1), p*q);
dataldot = zeros(N(1), p*q);
clear operators

operators = getMatrices(N(1), mu, "Burger");

disp("Start Burger Data Collection")

for i = 1:p*q

% High fidelity data
Xhtemp = semiEuler(operators.A, operators.H, x0(:, i), t);
Xhdottemp = getDot(Xhtemp, t);
Xhtemp = Xhtemp(:, 2:end);
datah(:, i) = Xhtemp;
datahdot(:, i) = Xhdottemp;

% Low fidelity data
Xltemp = semiEuler(operators.A, operators.H, x0l(:, i), tl);
Xltemp = interp1(tl, Xltemp', t)';
Xldottemp = getDot(Xltemp, t);
Xltemp = Xltemp(:, 2:end);
datal(:, i) = Xltemp;
dataldot(:, i) = Xldottemp;

% Tracking
if mod(i, p*q./100) == 0
    disp("Burger Percent Done: " + num2str(i*100./(p*q)) + "%")
end

end

save("BurgerEqData.mat","operators","datah", "datahdot", "datal", "dataldot", "para")
disp("Variables saved to BurgerEqData.mat");

%% Functions 

function param = getRandIC(p, q)
% INPUTS: 
%   p : number of random amplitude parameters
%   q : number of random frequency parameters
% OUTPUTS: 
%   param : (2 x pq) parameters [a; b] assoicated with each intial
%       condition
    if q ~= 0
    % Generating random amplitudes
    a = rand(p, 1);

    % Generating random frequencies
    b = randi(5, [q, 1]); %

    % parameters
    param = [repmat(a', 1, q); repelem(b', 1, p)];
    inds = randperm(p*q);
    param = param(:, inds);
    else 
        a = rand(p, 1);
        param = [a'; zeros(1, p)];
    end
end

function x0 = getIC(x, param)
% INPUTS: 
%   x : (N x 1) discretized spatial coordinate
%   param : (2 x pq) parameters [a; b] assoicated with each intial
%       condition
    x01 = param(1,:)'*sin(2.*pi.*x./x(end));
    x02 = cos(2*pi*param(2, :)'*x./x(end));
    x0 = (x01 + x02)';
end

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

function X = backwardEuler(A, x0, time)
% This function uses the backward euler's integration scheme.
% INPUTS:
%   A - (N x N) linear operator
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
    new_X = (eye(N) - dt*A)\(X(:, i));
    X(:, i+1) = new_X;
end
end

function X = semiEuler(A, H, x0, time)
% This function uses the semi-implicit euler's integration scheme.
% INPUTS:
%   A - (N x N) linear operator
%   H - (N x (N+1)C(2)) quadratic operator
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
    new_X = (eye(N) - dt.*A)\(X(:, i) + dt*H*kron(X(:, i), X(:, i)));
    X(:, i+1) = new_X;
end
end

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
end
end

