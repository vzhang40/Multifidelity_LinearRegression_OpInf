%% Problem Set-Up
clear; close all; clc

%% 1D Periodic Heat Equation
addpath("data_files")
addpath("functions")
load("heatEqData.mat")

% Getting POD basis with high fidelity state data
[U, ~, ~] = svd(datah, "econ");

% Budgets
p = [10, 50, 100, 1000, 2000];

% Costs, calculated using mesh dimension
w = [1; 16^2/128^2];

% Number of Replicates
reps = 100;

% Initializing data over basis sizes
storeAhatHF = cell(8, 1);
storeAhatMF = cell(8, 1);
storeerrorsHF = cell(8, 1);
storeerrorsMF = cell(8, 1);
opInfFull = cell(8, 1);

for r = 1:8
Ur = U(:, 1:r);

% Reduced Hifi Data
Xhr = Ur'*datah;
Xhrdot = Ur'*datahdot;

% Reduced Lofi Data
Xlr = Ur'*datal;
Xlrdot = Ur'*dataldot;

% Taking the covariance to be exact
XXT = (1./size(Xhr, 2))*(Xhr*Xhr'); 

% Full OpInf to be "True"
Ahat_true = (XXT+1e-3*eye(r))\((1./size(Xhr, 2))*Xhr*Xhrdot');
opInfFull{r} = Ahat_true;

%% Writing it in MFLR form

% Initializing Ahats
Ahat_hf = cell(length(p), reps);
Ahat_mf = cell(length(p), reps);

% Initializing Errors
hf_errors = zeros(length(p), reps);
mf_errors = zeros(length(p), reps);

% Finding control variate term
alpha = zeros(r, 1);
for i = 1:r
    ghigh = Xhr.*Xhrdot(i, :);
    glow = Xlr.*Xlrdot(i, :);

    ghighvar = ghigh - mean(ghigh);
    glowvar = glow - mean(glow);

    Gamhl = (ghighvar * glowvar') / (size(glowvar,2) - 1);
    Gamll = (glowvar * glowvar') / (size(glowvar,2) - 1);
    alpha(i) = trace(Gamhl)/trace(Gamll);
end

%% Find "optimal" sample size
rhos = zeros(r, 3);
m1 = zeros(r, length(p));
m2 = zeros(r, length(p));
rhos(:, 1) = ones(r, 1);
for i = 1:r
    temp = corrcoef(Xhrdot(i, :), Xlrdot(i, :));
    rhos(i, 1:end-1) = temp(1, :);
    rs = sqrt( w(1)*(rhos(i, 1:end-1).^2 - rhos(i, 2:end).^2) ./ (w'*(1 - rhos(i, 2).^2)))';
    m1(i, :) = floor(p./(w'*rs));
    m2(i, :) = floor(m1(i, :).*rs(2));
end
m2(m2 > size(datah, 2)) = size(datah, 2);
m = floor(p./w(1));

%% Looping for replicates
for k = 1:reps
    for j = 1:length(p)
    indsSF = randperm(size(datah, 2), m(j));
    indsMF = randperm(size(datah, 2), m2(r, j));

    % Initializing Ahats
    Ahat_hf{j, k} = zeros(r, r);
    Ahat_mf{j, k} = zeros(r, r);

    % Single fi Data
    X = Xhr(:, indsSF);
    Y = Xhrdot(:, indsSF)';

    % MFMC: Hifi Data 
    Xh = Xhr(:, indsMF(1:m1(r, j)));
    Yh = Xhrdot(:, indsMF(1:m1(r, j)))';

    % MFMC: Lofi Data
    Xl1 = Xlr(:, indsMF(1:m1(r, j)));
    Yl1 = Xlrdot(:, indsMF(1:m1(r, j)))';
    Xl2 = Xlr(:, indsMF);
    Yl2 = Xlrdot(:, indsMF)';
    
    for i = 1:r
        % single fidelity
        beta_single = (XXT+1e-3*eye(r))\((1./size(Y, 1))*X*Y(:, i));
        Ahat_hf{j, k}(i, :) = beta_single';

        % Multifidelity
        XYlm1 = (1./size(Yl1, 1))*Xl1*Yl1(:, i);
        XYlm2 = (1./size(Yl2, 1))*Xl2*Yl2(:, i);
        XYhm1 = (1./size(Yh, 1))*Xh*Yh(:, i);
        XYmulti = XYhm1 + alpha(i)*(XYlm2 - XYlm1);
        beta_multi = (XXT+1e-2*eye(r))\(XYmulti);
        Ahat_mf{j, k}(i, :) = beta_multi';
    end

    hf_errors(j, k) = norm(Ahat_hf{j, k} - Ahat_true, "fro")./norm(Ahat_true, "fro");
    mf_errors(j, k) = norm(Ahat_mf{j, k} - Ahat_true, "fro")./norm(Ahat_true, "fro");
    end

end 
    storeAhatHF{r} = Ahat_hf;
    storeAhatMF{r} = Ahat_mf;
    
    storeerrorsHF{r} = hf_errors;
    storeerrorsMF{r} = mf_errors;
end


%% Getting initial conditions
% Training Data
K = 5;
x = linspace(0, 1, size(U, 1));
x0 = getIC(x, para(:, 1));
t = linspace(0, 1 * K./1000, K + 1);

% getting true state
N = 128; mu = 0.01;
operators = getMatrices(N, mu, "heat");
A = operators.A;
true_state = backwardEuler(A, x0, t);

%
r = 8;
state_errorHF = cell(1, length(p));
state_errorMF = cell(1, length(p));
for i = 1:r
    Ur = U(:, 1:i);
    x0r = Ur'*x0;
    for j = 1:length(p)
        for k = 1:reps
            esti_rstateHF = backwardEuler(storeAhatHF{i}{j, k}, x0r, t);
            esti_rstateMF = backwardEuler(storeAhatMF{i}{j, k}, x0r, t);
            state_errorHF{j}(i, k) = norm(true_state - Ur*esti_rstateHF, "fro")./norm(true_state, "fro");
            state_errorMF{j}(i, k) = norm(true_state - Ur*esti_rstateMF, "fro")./norm(true_state, "fro");
        end
    end
end


% Plotting
rbg = orderedcolors("gem");
figure(4); clf(4)
set(gcf, 'Position', [100 100 800 450])
ranks = 1:6;
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
    stdHF = std(state_errorMF{r(i)}');
    stdMF = std(state_errorMF{r(i)}');
    stdHF = stdHF(ranks);
    stdMF = stdMF(ranks);
    patch([ranks'; flip(ranks)'], [meanHF-stdHF'; flip(meanHF+stdHF')], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
    patch([ranks'; flip(ranks)'], [meanMF-stdMF'; flip(meanMF+stdMF')], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
    xlabel("Basis Size")
    ylabel("State Error")
    sgtitle("Budget versus State Error: " + K + " Timesteps")
    title("Budget $p = " + p(r(i)) + "$")
    legend("Single Fidelity", "Multifidelity", "location", "southwest")
    xlim([min(ranks), max(ranks)])
    ylim([0.001, 1])
end


%% Plotting
rbg = orderedcolors("gem");
figure(1); clf(1)
for r = 1:8
    subplot(2, 4, r)
    meanHF = mean(storeerrorsHF{r}, 2);
    meanMF = mean(storeerrorsMF{r}, 2);
    loglog(p, meanHF, "Color", rbg(1, :))
    hold on
    loglog(p, meanMF, "Color", rbg(2, :))
    stdHF = std(storeerrorsMF{r}');
    stdMF = std(storeerrorsMF{r}');
    patch([p'; flip(p)'], [meanHF-stdHF'; flip(meanHF+stdHF')], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
    patch([p'; flip(p)'], [meanMF-stdMF'; flip(meanMF+stdMF')], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
    xlabel("Budget")
    ylabel("Error")
    sgtitle("Budget versus Operator Error")
    title("Basis Size $r = " + r + "$")
    legend("Single Fidelity", "Multifidelity")
    xlim([min(p), max(p)])
end


