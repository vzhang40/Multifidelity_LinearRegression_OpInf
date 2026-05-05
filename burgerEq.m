%% Problem Set-Up
clear; clc

%% 1D Periodic Heat Equation
addpath("data_files")
addpath("functions")
load("BurgerEqData.mat")

%% Analysis Set-up
% Budgets
p = [10, 50, 100, 1000, 2000];

% Costs, calculated using mesh dimension
w = [1; 16^2/128^2];

% Number of basis sizes tested
R = 12;

% Number of Replicates
reps = 100;

% Initializing data over basis sizes
storeOhatHF = cell(R, 1);
storeOhatMF = cell(R, 1);
storeerrorsHF_op = cell(R, 1);
storeerrorsMF_op = cell(R, 1);
storeerrorsHF_in = cell(R, 1);
storeerrorsMF_in = cell(R, 1);
storeerrorsMFh = cell(R, 1);
opInfs = cell(R, 1);
O_intrus = cell(R, 1);
alphas = cell(R,1);
ms = cell(R,1);
m1s = cell(R,1);
m2s = cell(R,1);
storeXXT = cell(R, 1);

% Getting POD data with hifi data
[U, S, ~] = svd(datah, "econ");

% Looping over basis sizes
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
storeXXT{r} = XXT;

% Full OpInf to be "True"
Ohat_opinf = opinf(Xhr_full, Xhrdot, size(Xhr, 2), XXT);
opInfs{r} = Ohat_opinf;

Ur_kron = zeros(size(Ur, 1).*(size(Ur, 1) + 1)./2, size(Ur,2));


% "Optimal" MFLR Parameters
% Finding control variate term
bigG = getStats(Xhr, Xhrdot, Xlr, Xlrdot);
alpha = optimalControlVariate(bigG);
[m1, m2, a1, a2] = mflrOptimalSampleSize(bigG, w, p); % multifidelity 
% [m1, m2] = mfmcOptimalSampleSize(p, w, Xhrdot, Xlrdot); % multifidelity

m2(m2 > size(datah, 2)) = size(datah, 2); % based on available data
% m2h(m2h > size(datah, 2)) = size(datah, 2); % based on available data

m = floor(p./w(1)); % single fidelity
% m1 = repmat(m, [r, 1]); m2 = m1 + 1;

% Initializing Ahats and operator errors
Ohat_hf = cell(length(p), reps);
Ohat_mf = cell(length(p), reps);
Ohat_mfh = cell(length(p), reps);
hf_errors_op = zeros(length(p), reps);
mf_errors_op = zeros(length(p), reps);
hf_errors_in = zeros(length(p), reps);
mf_errors_in = zeros(length(p), reps);
mf_errorsh = zeros(length(p), reps);

% Looping for replicates
for k = 1:reps
    for j = 1:length(p)
    % Initializing Ahats
    Ohat_hf{j, k} = opinf(Xhr_full, Xhrdot, m(j), XXT);
    Ohat_mf{j, k} =  mfopinf(Xhr_full, Xhrdot, Xlr_full, Xlrdot, m1(:, j), m2(:, j), alpha, XXT);

    hf_errors_op(j, k) = norm(Ohat_hf{j, k} - Ohat_opinf, "fro")./norm(Ohat_opinf, "fro");
    mf_errors_op(j, k) = norm(Ohat_mf{j, k} - Ohat_opinf, "fro")./norm(Ohat_opinf, "fro");
    end
    disp(k + " Rep out of 100")
end 
    storeOhatHF{r} = Ohat_hf;
    storeOhatMF{r} = Ohat_mf;
    
    storeerrorsHF_op{r} = hf_errors_op;
    storeerrorsMF_op{r} = mf_errors_op;
    disp(r + " Size Done")
end

% Plotting Operator Error
rbg = orderedcolors("gem");

figure; 
for r = 1:R
    subplot(2, ceil(R/2), r)
    loglog(p, mean(storeerrorsHF_op{r}, 2), "Color", rbg(1, :))
    hold on
    loglog(p, mean(storeerrorsMF_op{r}, 2), "Color", rbg(2, :))
    % loglog(p, mean(storeerrorsMFh_op{r}, 2), "Color", rbg(4, :))
    loglog(p, storeerrorsHF_op{r}, "Color", [rbg(1, :) 0.1], "LineWidth", 0.01)
    loglog(p, storeerrorsMF_op{r}, "Color", [rbg(2, :) 0.1], "LineWidth", 0.01)
    % loglog(p, storeerrorsMFh{r}, "Color", [rbg(4, :) 0.1], "LineWidth", 0.01)
    xlabel("Budget")
    ylabel("Error")
    sgtitle("Budget versus Operator Error Relative to OpInf")
    title("$r = " + r + "$")
    legend("Single Fidelity", "Multifidelity")
    xlim([min(p), max(p)])
end

%% Testing
t = linspace(0, 1, 1001);
% Full States

eHall = zeros(R, reps, length(p));
eMall = zeros(R, reps, length(p));
for k = 1:reps
    Xh = semiEuler(operators.A, operators.F, x0_test(:, k), t);
for i = 1:R
    Ur = U(:, 1:i);
    x0r = Ur'*x0_test(:, k);
    for j = 1:length(p)
        Oh = storeOhatHF{i}{j, k};
        Om = storeOhatMF{i}{j, k};
        esti_rstateHF = semiEuler(Oh(:, 1:i),Oh(:, i+1:end),  x0r, t);
        esti_rstateMF = semiEuler(Om(:, 1:i),Om(:, i+1:end),  x0r, t);
        
        sH = Ur*esti_rstateHF;
        sM = Ur*esti_rstateMF;
        eHall(i, k, j) = 0.001*norm(Xh - sH, "fro")./norm(Xh, "fro");
        eMall(i, k, j) = 0.001*norm(Xh - sM, "fro")./norm(Xh, "fro");
    end
end
disp(k + "% Done")
end

%
eHall(eHall >= 1) = NaN; eMall(eMall >= 1) = NaN;
eH = mean(eHall, 2, 'omitnan'); eH = reshape(eH, [R, length(p)]);
eM = mean(eMall, 2, 'omitnan'); eM = reshape(eM, [R, length(p)]);

%%
figure
subplot(1, 2, 1)
r = [1, 3, 5, 9, 12];

%
for i = 1:length(r)
loglog(p, eH(r(i), :), "-x", "Color", rbg(i, :))
hold on
end

for i = 1:length(r)
hold on
loglog(p, eM(r(i), :), "--o", "Color", rbg(i, :))
end
legendStrings = arrayfun(@(x) ['$r = ' num2str(x) '$'], r, 'UniformOutput', false);
legend(legendStrings, 'Interpreter', 'latex');
xlabel("Budget")
ylabel("Relative State Error")
title("State Error versus Budget")
%%
figure
subplot(1, 2, 1)
for i = 1:length(p)
semilogy(1:R, eH(:, i), "-x", "Color", rbg(i, :))
hold on
end
subplot(1, 2, 2)
for i = 1:length(p)
semilogy(1:R, eM(:, i), ":o", "Color", rbg(i, :))
hold on
end
xlabel("Dimension Size")
ylabel("Relative State Error")
title("Dimension Size versus Budget")
%% Getting initial conditions
% % Training Data

true_state = semiEuler(operators.A, operators.F, x0, t);
% 
% 
state_errorHF = cell(1, length(p));
state_errorMF = cell(1, length(p));

state_errorHF1 = cell(1, length(p));
state_errorMF1 = cell(1, length(p));
state_errorHF2 = cell(1, length(p));
state_errorMF2 = cell(1, length(p));
state_errorHF3 = cell(1, length(p));
state_errorMF3 = cell(1, length(p));

state_errorin1 = cell(1, length(p));
state_errorin2 = cell(1, length(p));
state_errorin3 = cell(1, length(p));
% 
% 

% Plotting State Errors
rbg = orderedcolors("gem");
figure(3); clf(3)
for r = 1:length(p)
    subplot(1, length(p), r)
    meanerrorHF1 = mean(state_errorHF1{r}, 3);
    meanerrorMF1 = mean(state_errorMF1{r}, 3);
    meanHF = mean(meanerrorHF1, 2);
    meanMF = mean(meanerrorMF1, 2);
    meanIN = mean(meanerrorIn1, 2);
    hold on
    semilogy(1:R, meanHF, "Color", rbg(1, :), "Marker", "o")
    semilogy(1:R, meanMF, "Color", rbg(2, :), "Marker", "+")
    xlabel("Basis Size")
    ylabel("Error")
    sgtitle("Budget versus State Error: Testing Data, Timestep 100 $t = 0.1$")
    title("Budget $p = " + p(r) + "$")
    legend("Single Fidelity", "Multifidelity")
    xlim([1, R])
end
%%
figure(4); clf(4)
for r = 1:length(p)
    subplot(1, length(p), r)
    meanerrorHF2 = mean(state_errorHF2{r}, 3);
    meanerrorMF2 = mean(state_errorMF2{r}, 3);
    meanHF = mean(meanerrorHF2, 2);
    meanMF = mean(meanerrorMF2, 2);
    hold on
    semilogy(1:R, meanHF, "Color", rbg(1, :), "Marker", "o")
    semilogy(1:R, meanMF, "Color", rbg(2, :), "Marker", "+")
    xlabel("Basis Size")
    ylabel("Error")
    sgtitle("Budget versus State Error: Testing Data, Timestep 200 $t = 0.2$")
    title("Budget $p = " + p(r) + "$")
    legend("Single Fidelity", "Multifidelity")
    xlim([1,R])
end
% 
%%
figure(5); clf(5)
for r = 1:length(p)
    subplot(1, length(p), r)
    meanerrorHF3 = mean(state_errorHF3{r}, 3);
    meanerrorMF3 = mean(state_errorMF3{r}, 3);
    meanHF = mean(meanerrorHF3, 2);
    meanMF = mean(meanerrorMF3, 2);
    hold on
    semilogy(1:R, meanHF, "Color", rbg(1, :), "Marker", "o")
    semilogy(1:R, meanMF, "Color", rbg(2, :), "Marker", "+")
    xlabel("Basis Size")
    ylabel("Error")
    sgtitle("Budget versus State Error: Testing Data, Timestep 500 ($t = 0.5$)")
    title("Budget $p = " + p(r) + "$")
    legend("Single Fidelity", "Multifidelity")
    xlim([1, R])
end
% % Plotting
rbg = orderedcolors("gem");
figure(4); clf(4)
set(gcf, 'Position', [100 100 800 450])
ranks = 1:R;
for i = 1:5
    r = 1:5;
    subplot(1, 5, i)
    meanHF = mean(state_errorHF{r(i)}, 2);
    meanMF = mean(state_errorMF{r(i)}, 2);
    meanHF = meanHF(ranks);
    meanMF = meanMF(ranks);
    semilogy(ranks, meanHF, "Color", rbg(1, :))
    hold on
    semilogy(ranks, meanMF, "Color", rbg(2, :))
    xlabel("Basis Size")
    ylabel("State Error")
    sgtitle("Budget versus State Error, Testing Data: " + K + " Timesteps")
    title("Budget $p = " + p(r(i)) + "$")
    legend("Single Fidelity", "Multifidelity", "location", "southwest")
    xlim([min(ranks), max(ranks)])
end
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



