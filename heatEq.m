%% Problem Set-Up
clear; close all; clc

%% 1D Periodic Heat Equation
addpath("data_files")
addpath("functions")
load("1heatEqData.mat")

%% Analysis Set-up
% Budgets
p = [10, 50, 100, 1000, 2000];

% Costs, calculated using mesh dimension
w = [1; 16^2/128^2];

% Number of basis sizes tested
R = 8;

% Number of Replicates
reps = 100;

% Initializing data over basis sizes
storeAhatHF = cell(R, 1);
storeAhatMF = cell(R, 1);
storeerrorsHF = cell(R, 1);
storeerrorsMF = cell(R, 1);
opInfFull = cell(R, 1);
alphas = cell(R,1);
ms = cell(R,1);
m1s = cell(R,1);
m2s = cell(R,1);

% Getting POD basis with high fidelity state data
[U, S, ~] = svd(datah, "econ");

%% Looping over basis sizes
for r = 1:R

% POD basis of size r
Ur = U(:, 1:r);

% Reducing high fidelity data
Xhr = Ur'*datah; Xhrdot = Ur'*datahdot;

% Reducing low fidelity data
Xlr = Ur'*datal; Xlrdot = Ur'*dataldot;

% Taking the covariance to be exact
XXT = (1./size(Xhr, 2))*(Xhr*Xhr'); 

% Full OpInf to be "True"
Ahat_true = opinf(Xhr, Xhrdot, size(Xhr, 2), XXT);
opInfFull{r} = Ahat_true;

%% "Optimal" MFLR Parameters
% Finding control variate term
bigG = getStats(Xhr, Xhrdot, Xlr, Xlrdot);
alpha = optimalControlVariate(bigG);
[m1, m2, a1, a2] = mflrOptimalSampleSize(bigG, w, p); % multifidelity 
% [m1, m2] = mfmcOptimalSampleSize(p, w, Xhrdot, Xlrdot); % multifidelity


%% Finding "optimal" sample size
m = floor(p./w(1)); % single fidelity

% Initializing Ahats and operator errors
Ahat_hf = cell(length(p), reps);
Ahat_mf = cell(length(p), reps);
hf_errors = zeros(length(p), reps);
mf_errors = zeros(length(p), reps);

%% Looping over replicates
for k = 1:reps
    for j = 1:length(p)
    Ahat_hf{j, k} = opinf(Xhr, Xhrdot, m(j), XXT);
    Ahat_mf{j, k} =  mfopinf(Xhr, Xhrdot, Xlr, Xlrdot, m1(:, j), m2(:, j), alpha, XXT);

    hf_errors(j, k) = norm(Ahat_hf{j, k} - Ahat_true, "fro")./norm(Ahat_true, "fro");
    mf_errors(j, k) = norm(Ahat_mf{j, k} - Ahat_true, "fro")./norm(Ahat_true, "fro");
    end
end 
% Storing variables
storeAhatHF{r} = Ahat_hf;
storeAhatMF{r} = Ahat_mf;
storeerrorsHF{r} = hf_errors;
storeerrorsMF{r} = mf_errors;
alphas{r} = alpha; ms{r} = m; m1s{r} = m1; m2s{r} = m2;
disp(r + "")
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
for r = 1:R
    subplot(2, ceil(R/2), r)
    hold on
    ylim([10^-3, 10^0])
end

%%
% figure(2); clf(2)
% rbg = orderedcolors("gem");
% for r = 1:R
%     meanHF = mean(storeerrorsHF{r}, 2);
%     meanMF = mean(storeerrorsMF{r}, 2);
% 
%     subplot(2, 3, r)
%     loglog(p, meanHF, "Color", rbg(1,:))
%     hold on
%     loglog(p, meanMF,"Color", rbg(2,:))
%     loglog(p, storeerrorsHF{r}, ...
%         "Color", [rbg(1,:) 0.1], "LineWidth", 0.1)
%     loglog(p, storeerrorsMF{r}, ...
%         "Color", [rbg(2,:) 0.1], "LineWidth", 0.1) 
%     xlabel("Budget")
%     ylabel("Error")
%     title("Basis Size $r = " + r + "$", "Interpreter", "latex")
%     xlim([min(p), max(p)])
%     legend("Single Fidelity", "Multi-Fidelity ")
% end
% sgtitle("Budget versus Operator Error")
% 
% 
% %% Getting initial conditions
% % Training Data
% K = 100; % max timesteps to test
% numIC = 100;
% x = linspace(0, 1, size(U, 1));
% x0 = getIC(x, para(:, 1:numIC));
% t = linspace(0, 1 * K./1000, K + 1);
% 
% % getting true state
% N = 128; mu = 0.01;
% operators = getMatrices(N, mu, "heat");
% A = operators.A;
% for m = 1:numIC
% true_state = backwardEuler(A, x0(:, m), t);
% state_errorHF1 = cell(1, length(p));
% state_errorMF1 = cell(1, length(p));
% state_errorHF2 = cell(1, length(p));
% state_errorMF2 = cell(1, length(p));
% state_errorHF3 = cell(1, length(p));
% state_errorMF3 = cell(1, length(p));
% for j = 1:R
%     Ur = U(:, 1:j);
%     x0r = Ur'*x0(:, m);
%     for i = 1:length(p)
%         for k = 1:reps
%             esti_rstateHF = backwardEuler(storeAhatHF{j}{i, k}, x0r, t);
%             esti_rstateMF = backwardEuler(storeAhatMF{j}{i, k}, x0r, t);
% 
%             state_errorHF1{i}(j, m, k) = norm(true_state(:, 5) - Ur*esti_rstateHF(:, 5), "fro")./norm(true_state(:, 5), "fro");
%             state_errorMF1{i}(j, m, k) = norm(true_state(:, 5) - Ur*esti_rstateMF(:, 5), "fro")./norm(true_state(:, 5), "fro");
% 
%             state_errorHF2{i}(j, m, k) = norm(true_state(:, 20) - Ur*esti_rstateHF(:, 20), "fro")./norm(true_state(:, 20), "fro");
%             state_errorMF2{i}(j, m, k) = norm(true_state(:, 20) - Ur*esti_rstateMF(:, 20), "fro")./norm(true_state(:, 20), "fro");
% 
%             state_errorHF3{i}(j, m, k) = norm(true_state(:, 100) - Ur*esti_rstateHF(:, 100), "fro")./norm(true_state(:, 100), "fro");
%             state_errorMF3{i}(j, m, k) = norm(true_state(:, 100) - Ur*esti_rstateMF(:, 100), "fro")./norm(true_state(:, 100), "fro");
%         end
%     end
% end
% end  
% %%
% % Plotting State Errors
% rbg = orderedcolors("gem");
% figure(3); clf(3)
% for r = 1:length(p)
%     subplot(1, length(p), r)
%     meanerrorHF1 = mean(state_errorHF1{r}, 3);
%     meanerrorMF1 = mean(state_errorMF1{r}, 3);
%     meanHF = mean(meanerrorHF1, 2);
%     meanMF = mean(meanerrorMF1, 2);
%     semilogy(1:R, meanHF, "Color", rbg(1, :))
%     hold on
%     semilogy(1:R, meanMF, "Color", rbg(2, :))
%     for i = 1:100
%     semilogy(1:R, meanerrorHF1(:, i), ...
%         "Color", [rbg(1,:) 0.2], "LineWidth", 0.1)
%     semilogy(1:R, meanerrorMF1(:, i), ...
%         "Color", [rbg(2,:) 0.2], "LineWidth", 0.1) 
%     end
%     xlabel("Basis Size")
%     ylabel("Error")
%     sgtitle("Budget versus State Error: Training Data, Timestep 5")
%     title("Budget $p = " + p(r) + "$")
%     legend("Single Fidelity", "Multifidelity")
%     xlim([1, 6])
% end
% %%
% figure(4); clf(4)
% for r = 1:length(p)
%     subplot(1, length(p), r)
%     meanerrorHF2 = mean(state_errorHF2{r}, 3);
%     meanerrorMF2 = mean(state_errorMF2{r}, 3);
%     meanHF = mean(meanerrorHF2, 2);
%     meanMF = mean(meanerrorMF2, 2);
%     semilogy(1:R, meanHF, "Color", rbg(1, :))
%     hold on
%     semilogy(1:R, meanMF, "Color", rbg(2, :))
%     loglog(1:R, meanerrorHF2, ...
%         "Color", [rbg(1,:) 0.1], "LineWidth", 0.1)
%     loglog(1:R, meanerrorMF2, ...
%         "Color", [rbg(2,:) 0.1], "LineWidth", 0.1) 
%     xlabel("Basis Size")
%     ylabel("Error")
%     sgtitle("Budget versus State Error: Training Data, Timestep 20")
%     title("Budget $p = " + p(r) + "$")
%     legend("Single Fidelity", "Multifidelity")
%     xlim([1, 6])
%     ylim([10^-5, 10^-1])
% end
% 
% %%
% figure(5); clf(5)
% for r = 1:length(p)
%     subplot(1, length(p), r)
%     meanerrorHF3 = mean(state_errorHF3{r}, 3);
%     meanerrorMF3 = mean(state_errorMF3{r}, 3);
%     meanHF = mean(meanerrorHF3, 2);
%     meanMF = mean(meanerrorMF3, 2);
%     semilogy(1:R, meanHF, "Color", rbg(1, :))
%     hold on
%     semilogy(1:R, meanMF, "Color", rbg(2, :))
%     loglog(1:R, meanerrorHF2, ...
%         "Color", [rbg(1,:) 0.1], "LineWidth", 0.1)
%     loglog(1:R, meanerrorMF2, ...
%         "Color", [rbg(2,:) 0.1], "LineWidth", 0.1) 
%     xlabel("Basis Size")
%     ylabel("Error")
%     sgtitle("Budget versus State Error: Training Data, Timestep 100")
%     title("Budget $p = " + p(r) + "$")
%     legend("Single Fidelity", "Multifidelity")
%     xlim([1, 6])
%     ylim([10^-5, 10^-1])
% end
% 
% 
% 
% 
% 





% % 
% % %%
% % rbg = orderedcolors("gem");
% % figure(3); clf(3)
% % for r = 1:length(p)
% %     subplot(1, length(p), r)
% %     meanerrorHF2 = mean(state_errorHF2{r}, 3);
% %     meanerrorMF2 = mean(state_errorMF2{r}, 3);
% %     meanHF = mean(meanerrorHF2, 2);
% %     meanMF = mean(meanerrorMF2, 2);
% %     semilogy(1:R, meanHF, "Color", rbg(1, :))
% %     hold on
% %     semilogy(1:R, meanMF, "Color", rbg(2, :))
% %     stdHF = std(meanerrorHF2');
% %     stdMF = std(meanerrorHF2');
% %     patch([(1:R)'; (1:R)'], [meanHF-stdHF'; flip(meanHF+stdHF')], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
% %     patch([(1:R)'; (1:R)'], [meanMF-stdMF'; flip(meanMF+stdMF')], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
% %     xlabel("Basis Size")
% %     ylabel("Error")
% %     sgtitle("Budget versus State Error: Training Data, Timestep 20")
% %     title("Budget $p = " + p(r) + "$")
% %     legend("Single Fidelity", "Multifidelity")
% %     xlim([1, 6])
% %     ylim([10^-5, 10^-1])
% % end
% % 
% % %%
% % rbg = orderedcolors("gem");
% % figure(4); clf(4)
% % for r = 1:length(p)
% %     subplot(1, length(p), r)
% %     meanerrorHF3 = mean(state_errorHF3{r}, 3);
% %     meanerrorMF3 = mean(state_errorMF3{r}, 3);
% %     meanHF = mean(meanerrorHF3, 2);
% %     meanMF = mean(meanerrorMF3, 2);
% %     semilogy(1:R, meanHF, "Color", rbg(1, :))
% %     hold on
% %     semilogy(1:R, meanMF, "Color", rbg(2, :))
% %     stdHF = std(meanerrorHF3');
% %     stdMF = std(meanerrorHF3');
% %     patch([(1:R)'; (1:R)'], [meanHF-stdHF'; flip(meanHF+stdHF')], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
% %     patch([(1:R)'; (1:R)'], [meanMF-stdMF'; flip(meanMF+stdMF')], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
% %     xlabel("Basis Size")
% %     ylabel("Error")
% %     sgtitle("Budget versus State Error: Training Data, Timestep 100")
% %     title("Budget $p = " + p(r) + "$")
% %     legend("Single Fidelity", "Multifidelity")
% %     xlim([1, 6])
% %     ylim([10^-5, 10^-1])
% % end

% for r = 1:6
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
