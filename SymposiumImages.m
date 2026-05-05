clear; clc
addpath("functions\")

%% Model Dimensions and Parameters
L = 1; % length
T = 5.5; % final time
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



num_samples = 1;
max_param = 60;

figure(67);
for i = 1:3
xi = getRandIC(max_param, num_samples);
x0 = getIC(x, xi);
x0l = getIC(xl, xi);

operators = getMatrices(N(1), mu, "Burger");
operatorsl = getMatrices(N(2), mu, "Burger");

% High fidelity data
Xh = semiEuler(operators.A, operators.F, x0, t);

% Low fidelity data
Xl = semiEuler(operatorsl.A, operatorsl.F, x0l, tl);


[Th, Xh_grid] = meshgrid(x, t);
[Tl, Xl_grid] = meshgrid(xl, tl);


%% ---- Fine Mesh Heatmap ----
subplot(3,1,i)
imagesc(t, x, Xh)
set(gca, 'XTickLabel', [], 'YTickLabel', [])
colormap parula
colorbar
if i == 1
title("1D Viscous Burger's Equation Solution")
end
if i == 3
xlabel('$t$')
end
ylabel('$\eta$')
axis tight
axis equal
end
% 
% %% ---- Coarse Mesh Heatmap ----
% subplot(1,2,2)
% imagesc(tl, xl, Xl)
% set(gca, 'XTickLabel', [], 'YTickLabel', [])
% colormap parula
% colorbar
% title('Low Fidelity PDE')
% xlabel('$t$'), ylabel('$\eta$')
% axis tight
% axis equal

%%