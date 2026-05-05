% This script collects 120000 snapshot vectors for the 1D Heat Equation and
% the 1D Viscous Burger's Equation for a bifidelity operator inference
% experiment.

clear; close all; clc;
addpath("functions\")

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

num_samples = 200000;
max_param = 60;
xi = getRandIC(max_param, num_samples);
x0 = getIC(x, xi);
x0l = getIC(xl, xi);
x0l = interp1(xl, x0l, x, "linear"); 

num_test = 100;
xi_test = getRandIC(max_param, num_test);
x0_test = getIC(x, xi_test);

datah = zeros(N(1), num_samples);
datahdot = zeros(N(1), num_samples);
datal = zeros(N(1), num_samples);
dataldot = zeros(N(1), num_samples);

%% Generating HEAT EQ Data
operators = getMatrices(N(1), mu, "heat");

disp("Start Heat Data Collection")
for i = 1:num_samples

% High fidelity data
Xhtemp = backwardEuler(operators.A, x0(:, i), t);
Xhdottemp = getDot(Xhtemp, t);
Xhtemp = Xhtemp(:, 1:end-1);
datah(:, i) = Xhtemp;
datahdot(:, i) = Xhdottemp;

% Low fidelity data
Xltemp = backwardEuler(operators.A, x0l(:, i), tl);
Xltemp = interp1(tl, Xltemp', t)';
Xldottemp = getDot(Xltemp, t);
Xltemp = Xltemp(:, 1:end-1);
datal(:, i) = Xltemp;
dataldot(:, i) = Xldottemp;

% Tracking
if mod(i, num_samples./20) == 0
    disp("Heat Percent Done: " + num2str(i*100./num_samples) + "%")
end

end

save("heatEqData.mat","operators","datah", "datahdot", "datal", "dataldot", "x0", "x0_test")
disp("Variables saved to new_heatEqData.mat")

%% Generating BURGER EQ Data

% initializing
datah = zeros(N(1), num_samples);
datahdot = zeros(N(1), num_samples);
datal = zeros(N(1), num_samples);
dataldot = zeros(N(1), num_samples);
clear operators

operators = getMatrices(N(1), mu, "Burger");

disp("Start Burger Data Collection")

for i = 1:num_samples

% High fidelity data
Xhtemp = semiEuler(operators.A, operators.F, x0(:, i), t);
Xhdottemp = getDot(Xhtemp, t);
Xhtemp = Xhtemp(:, 1:end-1);
datah(:, i) = Xhtemp;
datahdot(:, i) = Xhdottemp;

% Low fidelity data
Xltemp = semiEuler(operators.A, operators.F, x0l(:, i), tl);
Xltemp = interp1(tl, Xltemp', t)';
Xldottemp = getDot(Xltemp, t);
Xltemp = Xltemp(:, 1:end-1);
datal(:, i) = Xltemp;
dataldot(:, i) = Xldottemp;

% Tracking
if mod(i, num_samples./100) == 0
    disp("Burger Percent Done: " + num2str(i*100./ num_samples) + "%")
end

end

save("new_BurgerEqData.mat","operators","datah", "datahdot", "datal", "dataldot", "x0", "x0_test")
disp("Variables saved to new_BurgerEqData.mat");



