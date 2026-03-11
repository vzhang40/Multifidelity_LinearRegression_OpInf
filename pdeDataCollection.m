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
if mod(i, p*q./20) == 0
    disp("Heat Percent Done: " + num2str(i*100./(p*q)) + "%")
end

end

save("1heatEqData.mat","operators","datah", "datahdot", "datal", "dataldot", "para")
disp("Variables saved to 1heatEqData.mat")

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
if mod(i, p*q./100) == 0
    disp("Burger Percent Done: " + num2str(i*100./(p*q)) + "%")
end

end

save("1BurgerEqData.mat","operators","datah", "datahdot", "datal", "dataldot", "para")
disp("Variables saved to 1BurgerEqData.mat");



