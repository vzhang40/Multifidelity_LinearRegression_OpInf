%% Problem Set-Up
clear; close all; clc

%% 1D Periodic Heat Equation
load("heatEqData.mat")

% Getting POD data with hifi data
[U, ~, ~] = svd(datah, "econ");

% Budgets
p = [10, 50, 100, 1000, 5000];

% Costs
w1 = 1; w2 = 0.01;

% Single fidelity sample sizes
m = p./w1;

% Multifidelity sample sizes
m1 = 0.8*p;
m2 = (p - m1)./w2;

% Number of Replicates
reps = 100;

storerhos = cell(8, 1);
sigma1s = cell(8, 1);
sigma2s = cell(8, 1);
alphas = cell(8, 1);

storeAhatHF = cell(8, 1);
storeAhatMF = cell(8, 1);
storeAhatMF1 = cell(8, 1);
storeAhatMF2 = cell(8, 1);

storeerrorsHF = cell(8, 1);
storeerrorsMF = cell(8, 1);
storeerrorsMF1 = cell(8, 1);
storeerrorsMF2 = cell(8, 1);
%%
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

%% Writing it in MFLR form

% Initializing Ahats
Ahat_hf = cell(length(p), reps);
Ahat_mf = cell(length(p), reps);
Ahat_mf1 = cell(length(p), reps);
Ahat_mf2 = cell(length(p), reps);

% Initializing Errors
hf_errors = zeros(length(p), reps);
mf_errors = zeros(length(p), reps);
mf_errors1 = zeros(length(p), reps); 
mf_errors2 = zeros(length(p), reps); 

% % 
rhos = zeros(r, 1);
sigma1 = zeros(r, 1);
sigma2 = zeros(r, 1);
for i = 1:r
    rho = corrcoef(Xhrdot(i, :)', Xlrdot(i, :)');
    rhos(i)= rho(1, 2);
    sigma1(i) = std(Xhrdot(i, :)');
    sigma2(i) = std(Xlrdot(i, :)');
end
alpha1 = rhos.*sigma1./sigma2;

alpha = zeros(r, 1);
alpha2 = cell(r, 1);
for i = 1:r
    ghigh = Xhr.*Xhrdot(i, :);
    glow = Xlr.*Xlrdot(i, :);

    ghighvar = ghigh - mean(ghigh);
    glowvar = glow - mean(glow);

    Gamhl = (ghighvar * glowvar') / (size(glowvar,2) - 1);
    Gamll =(glowvar * glowvar') / (size(glowvar,2) - 1);
    alpha(i) = trace(Gamhl)/trace(Gamll);
    alpha2{i} = (Gamhl + 1e-3*eye(r))\Gamll;
end

% storerhos{r} = rhos;
% sigma1s{r} = sigma1;
% sigma2s{r} = sigma2;
% alphas{r} = alpha;
    
for k = 1:reps

    for j = 1:length(p)
    indsSF = randperm(size(datah, 2), m(j));
    indsMF = randperm(size(datah, 2), m2(j));

    % Initializing Ahats
    Ahat_hf{j, k} = zeros(r, r);
    Ahat_mf{j, k} = zeros(r, r);

    % Single fi Data
    X = Xhr(:, indsSF);
    Y = Xhrdot(:, indsSF)';

    % MFMC: Hifi Data 
    Xh = Xhr(:, indsMF(1:m1(j)));
    Yh = Xhrdot(:, indsMF(1:m1(j)))';

    % MFMC: Lofi Data
    Xl1 = Xlr(:, indsMF(1:m1(j)));
    Yl1 = Xlrdot(:, indsMF(1:m1(j)))';
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

        XYmulti1 = XYhm1 + alpha1(i)*(XYlm2 - XYlm1);
        beta_multi1 = (XXT+1e-2*eye(r))\(XYmulti1);
        Ahat_mf1{j, k}(i, :) = beta_multi1';

        XYmulti2 = XYhm1 + alpha2{i}*(XYlm2 - XYlm1);
        beta_multi2 = (XXT+1e-2*eye(r))\(XYmulti2);
        Ahat_mf2{j, k}(i, :) = beta_multi2';
    end

    hf_errors(j, k) = norm(Ahat_hf{j, k} - Ahat_true, "fro");%./norm(Ahat_true, "fro");
    mf_errors(j, k) = norm(Ahat_mf{j, k} - Ahat_true, "fro");%./norm(Ahat_true, "fro");
    mf_errors1(j, k) = norm(Ahat_mf1{j, k} - Ahat_true, "fro");%./norm(Ahat_true, "fro");
    mf_errors2(j, k) = norm(Ahat_mf2{j, k} - Ahat_true, "fro");%./norm(Ahat_true, "fro");
    end

end 

    storeAhatMF{r} = Ahat_hf;
    storeAhatMF{r} = Ahat_mf;
    storeAhatMF1{r} = Ahat_mf1;
    storeAhatMF2{r} = Ahat_mf2;
    
    storeerrorsHF{r} = hf_errors;
    storeerrorsMF{r} = mf_errors;
    storeerrorsMF1{r} = mf_errors1;
    storeerrorsMF2{r} = mf_errors2;
end

%% Plotting
rbg = orderedcolors("gem");
figure(1)
for r = 1:8
    subplot(2, 4, r)
    loglog(p, mean(storeerrorsHF{r}, 2), "Color", rbg(1, :))
    hold on
    loglog(p, mean(storeerrorsMF{r}, 2), "Color", rbg(2, :))
    loglog(p, storeerrorsHF{r}, "Color", [rbg(1, :) 0.1], "LineWidth", 0.01)
    loglog(p, storeerrorsMF{r}, "Color", [rbg(2, :) 0.1], "LineWidth", 0.01)
    xlabel("Budget")
    ylabel("Error")
    sgtitle("Budget versus Operator Error Relative to Optimal OpInf")
    title("Basis Size $r = " + r + "$")
    legend("Single Fidelity", "Multifidelity")
    xlim([min(p), max(p)])
end

%%
figure(2); clf(2)
for r = 1:8
    subplot(2, 4, r)
    meanHF = mean(storeerrorsHF{r}, 2);
    meanMF = mean(storeerrorsMF{r}, 2);
    meanMF1 = mean(storeerrorsMF1{r}, 2);
    meanMF2 = mean(storeerrorsMF2{r}, 2);
    loglog(p, meanHF, "Color", rbg(1, :))
    hold on
    loglog(p, meanMF, "Color", rbg(2, :))
    loglog(p, meanMF1, "Color", rbg(3, :))
    loglog(p, meanMF2, "Color", rbg(4, :))
    stdHF = std(storeerrorsMF{r}');
    stdMF = std(storeerrorsMF{r}');
    stdMF1 = std(storeerrorsMF1{r}');
    stdMF2 = std(storeerrorsMF2{r}');
    % patch([p'; flip(p)'], [meanHF-stdHF'; flip(meanHF+stdHF')], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
    % patch([p'; flip(p)'], [meanMF-stdMF'; flip(meanMF+stdMF')], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
    % patch([p'; flip(p)'], [meanMF1-stdMF1'; flip(meanMF1+stdMF1')], 'y', 'FaceAlpha', 0.1, 'EdgeColor','none')
    % patch([p'; flip(p)'], [meanMF2-stdMF2'; flip(meanMF2+stdMF2')], 'p', 'FaceAlpha', 0.1, 'EdgeColor','none')
    xlabel("Budget")
    ylabel("Error")
    sgtitle("Budget versus Operator Error")
    title("Basis Size $r = " + r + "$")
    legend("Single Fidelity", "Multifidelity: Optimal Scalar", "Multifidelity: Heuristic", "Multifidelity: Optimal Matrix")
    xlim([min(p), max(p)])
end

%%
figure(3); clf(3)
for r = 1:8
    subplot(2, 4, r)
    meanHF = mean(storeerrorsHF{r}, 2);
    meanMF = mean(storeerrorsMF{r}, 2);
    loglog(p, meanHF, "Color", rbg(1, :))
    hold on
    loglog(p, meanMF, "Color", rbg(2, :))
    stdHF = std(storeerrorsMF{r}');
    stdMF = std(storeerrorsMF{r}');
    stdMF1 = std(storeerrorsMF1{r}');
    patch([p'; flip(p)'], [meanHF-stdHF'; flip(meanHF+stdHF')], 'b', 'FaceAlpha', 0.1, 'EdgeColor','none')
    patch([p'; flip(p)'], [meanMF-stdMF'; flip(meanMF+stdMF')], 'r', 'FaceAlpha', 0.1, 'EdgeColor','none')
    xlabel("Budget")
    ylabel("Error")
    sgtitle("Budget versus Operator Error")
    title("Basis Size $r = " + r + "$")
    legend("Single Fidelity", "Multifidelity")
    xlim([min(p), max(p)])
end

