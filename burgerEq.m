%% Problem Set-Up
clear; close all; clc

%% 1D Periodic Heat Equation
addpath("data_files")
addpath("functions")
load("1BurgerEqData.mat")

%% Analysis Set-up
% Budgets
p = [10, 50, 100, 1000, 2000, 3000];

% Costs, calculated using mesh dimension
w = [1; 16^2/128^2];

% Number of basis sizes tested
R = 6;

% Number of Replicates
reps = 100;

% Initializing data over basis sizes
storeOhatHF = cell(R, 1);
storeOhatMF = cell(R, 1);
storeerrorsHF = cell(R, 1);
storeerrorsMF = cell(R, 1);
opInfFull = cell(R, 1);
alphas = cell(R,1);
ms = cell(R,1);
m1s = cell(R,1);
m2s = cell(R,1);

% Getting POD data with hifi data
[U, S, ~] = svd(datah, "econ");

%% Looping over basis sizes
for r = 1:R
Ur = U(:, 1:r);

% Reduced Hifi Data
Xhr = Ur'*datah;
Xhrdot = Ur'*datahdot;

Xhr_kron = zeros(size(Xhr, 1).*(size(Xhr, 1) + 1)./2, size(Xhr,2));
for j = 1:size(Xhr, 2)
    Xhr_kron(:, j) = uniKron(Xhr(:, j));
end

Xhr_full = [Xhr' Xhr_kron']';

% Reduced Lofi Data
Xlr = Ur'*datal;
Xlrdot = Ur'*dataldot;

Xlr_kron = zeros(size(Xlr, 1).*(size(Xlr, 1) + 1)./2, size(Xlr,2));
for j = 1:size(Xlr, 2)
    Xlr_kron(:, j) = uniKron(Xlr(:, j));
end

Xlr_full = [Xlr' Xlr_kron']';

% Taking the covariance to be exact
XXT = (1./size(Xhr, 2))*(Xhr_full*Xhr_full'); 

% Full OpInf to be "True"
Ohat_true = opinf(Xhr_full, Xhrdot, size(Xhr, 2), XXT);
opInfFull{r} = Ohat_true;

%% "Optimal" MFLR Parameters
% Finding control variate term
bigG = getStats(Xhr, Xhrdot, Xlr, Xlrdot);
alpha = optimalControlVariate(bigG);
%[m1, m2, a1, a2] = mflrOptimalSampleSize(bigG, w, p); % multifidelity 
% [m1, m2] = mfmcOptimalSampleSize(p, w, Xhrdot, Xlrdot); % multifidelity
% m2(m2 > size(datah, 2)) = size(datah, 2); % based on available data
m = floor(p./w(1)); % single fidelity
m1 = repmat(m, [r, 1]);
m2 = m1;

% Initializing Ahats and operator errors
Ohat_hf = cell(length(p), reps);
Ohat_mf = cell(length(p), reps);
hf_errors = zeros(length(p), reps);
mf_errors = zeros(length(p), reps);

%% Looping for replicates
for k = 1:reps
    for j = 1:length(p)
    % Initializing Ahats
    Ohat_hf{j, k} = opinf(Xhr_full, Xhrdot, m(j), XXT);
    Ohat_mf{j, k} =  mfopinf(Xhr_full, Xhrdot, Xlr_full, Xlrdot, m1(:, j), m2(:, j), alpha, XXT);

    hf_errors(j, k) = norm(Ohat_hf{j, k} - Ohat_true, "fro")./norm(Ohat_true, "fro");
    mf_errors(j, k) = norm(Ohat_mf{j, k} - Ohat_true, "fro")./norm(Ohat_true, "fro");
    disp(k)
    end
end 
    storeOhatHF{r} = Ohat_hf;
    storeOhatMF{r} = Ohat_mf;
    
    storeerrorsHF{r} = hf_errors;
    storeerrorsMF{r} = mf_errors;
end

%% Plotting Operator Error
rbg = orderedcolors("gem");

figure(2); clf(2)
for r = 1:R
    subplot(2, ceil(R/2), r)
    loglog(p, mean(storeerrorsHF{r}, 2), "Color", rbg(1, :))
    hold on
    loglog(p, mean(storeerrorsMF{r}, 2), "Color", rbg(2, :))
    loglog(p, storeerrorsHF{r}, "Color", [rbg(1, :) 0.1], "LineWidth", 0.01)
    loglog(p, storeerrorsMF{r}, "Color", [rbg(2, :) 0.1], "LineWidth", 0.01)
    xlabel("Budget")
    ylabel("Error")
    % sgtitle("Budget versus Operator Error Relative to Optimal OpInf")
    title("Reduced Dimension $r = " + r + "$")
    legend("Single Fidelity", "Multifidelity")
    xlim([min(p), max(p)])
end

%%
% %% Getting initial conditions
% % Training Data
% K = 1000;
% x = linspace(0, 1, size(U, 1));
% x0 = getIC(x, para(:, 1));
% t = linspace(0, 1 * K./1000, K + 1);
% 
% % getting true state
% N = 128; mu = 0.01;
% operators = getMatrices(N, mu, "Burger");
% A = operators.A;
% F = operators.F;
% true_state = semiEuler(A, F, x0, t);
% 
% 
% state_errorHF = cell(1, length(p));
% state_errorMF = cell(1, length(p));
% for i = 1:2
%     Ur = U(:, 1:i);
%     x0r = Ur'*x0;
%     for j = 1:length(p)
%         for k = 1:reps
%             esti_rstateHF = semiEuler(storeOhatHF{i}{j, k}(:, 1:i),storeOhatHF{i}{j, k}(:, i+1:end),  x0r, t);
%             esti_rstateMF = semiEuler(storeOhatMF{i}{j, k}(:, 1:i),storeOhatMF{i}{j, k}(:, i+1:end), x0r, t);
%             state_errorHF{j}(i, k) = norm(true_state - Ur*esti_rstateHF, "fro")./norm(true_state, "fro");
%             state_errorMF{j}(i, k) = norm(true_state - Ur*esti_rstateMF, "fro")./norm(true_state, "fro");
%         end
%     end
% end
% 
% 
% % Plotting
% rbg = orderedcolors("gem");
% figure(4); clf(4)
% set(gcf, 'Position', [100 100 800 450])
% ranks = 1:2;
% for i = 1:5
%     r = 1:5;
%     subplot(1, 5, i)
%     meanHF = mean(state_errorHF{r(i)}, 2);
%     meanMF = mean(state_errorMF{r(i)}, 2);
%     meanHF = meanHF(ranks);
%     meanMF = meanMF(ranks);
%     semilogy(ranks, meanHF, "Color", rbg(1, :))
%     hold on
%     semilogy(ranks, meanMF, "Color", rbg(2, :))
%     stdHF = std(state_errorMF{r(i)}');
%     stdMF = std(state_errorMF{r(i)}');
%     stdHF = stdHF(ranks);
%     stdMF = stdMF(ranks);
%     patch([ranks'; flip(ranks)'], [meanHF-stdHF'; flip(meanHF+stdHF')], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
%     patch([ranks'; flip(ranks)'], [meanMF-stdMF'; flip(meanMF+stdMF')], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
%     xlabel("Basis Size")
%     ylabel("State Error")
%     sgtitle("Budget versus State Error: " + K + " Timesteps")
%     title("Budget $p = " + p(r(i)) + "$")
%     legend("Single Fidelity", "Multifidelity", "location", "southwest")
%     xlim([min(ranks), max(ranks)])
% end
% 
% 
% %% Plotting
% rbg = orderedcolors("gem");
% figure(2); clf(2)
% 
% for r = 1:R
%     meanHF = mean(storeerrorsHF{r}, 2);
%     meanMF = mean(storeerrorsMF{r}, 2);
% 
%     subplot(2, 3, r)
%     loglog(p, meanHF, ...
%         "Color", rbg(1,:))
%     hold on
% 
%     loglog(p, meanMF, ...
%         "Color", rbg(2,:))
% 
%     % --- Plot all Single-Fidelity trajectories (faint) ---
%     loglog(p, storeerrorsHF{r}, ...
%         "Color", [rbg(1,:) 0.1], "LineWidth", 0.1)   % transparent
% 
%     % --- Plot all Multi-Fidelity trajectories (faint) ---
%     loglog(p, storeerrorsMF{r}, ...
%         "Color", [rbg(2,:) 0.1], "LineWidth", 0.1)   % transparent
% 
%     xlabel("Budget")
%     ylabel("Error")
%     title("Basis Size $r = " + r + "$", "Interpreter", "latex")
%     xlim([min(p), max(p)])
%     ylim([10^-5, 10^4])
%     legend("Single Fidelity", "Multi-Fidelity", "Location", "south")
%     axis square
% end
% 
% sgtitle("Budget versus Operator Error")



% for r = 1:R
%     subplot(2, 3, r)
%     meanHF = mean(storeerrorsHF{r}, 2);
%     meanMF = mean(storeerrorsMF{r}, 2);
%     loglog(p, meanHF, "Color", rbg(1, :))
%     hold on
%     loglog(p, meanMF, "Color", rbg(2, :))
%     stdHF = std(storeerrorsMF{r}');
%     stdMF = std(storeerrorsMF{r}');
%     patch([p'; flip(p)'], [meanHF-stdHF'; flip(meanHF+stdHF')], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
%     patch([p'; flip(p)'], [meanMF-stdMF'; flip(meanMF+stdMF')], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
%     xlabel("Budget")
%     ylabel("Error")
%     sgtitle("Budget versus Operator Error")
%     title("Basis Size $r = " + r + "$")
%     legend("Single Fidelity", "Multifidelity")
%     xlim([min(p), max(p)])
% end



