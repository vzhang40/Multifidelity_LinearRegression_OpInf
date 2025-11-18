%% Problem Set-Up
clear; close all; clc

%% Data
load("data1.mat")

% Getting POD data with hifi data
[U, ~, ~] = svd(datah, "econ");

% Basis Size
figure(1)
for r = 2:9
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
Ahat_true = Xhr'\Xhrdot';

%% Writing it in MFLR form
% Testing Basis Size 8
p = [10, 50, 100, 1000];

% Single fidelity sample sizes
m = p;

% Multifidelity sample sizes
m1 = [8, 40, 80, 800];
m2 = (p - m1).*100;

% Initializing Ahats
Ahat_hf = cell(3, 1);
Ahat_mf = cell(3, 1);

% Initializing Errors
hf_errors = zeros(3, 1);
mf_errors = zeros(3, 1);


rhos = zeros(r, 1);
for i = 1:r
    rho = corrcoef(Xhr*Xhrdot(i, :)', Xlr(:, 1:500)*Xlrdot(i, 1:500)');
    rhos(i) = rho(1, 2);
end

for j = 1:length(p)

    % Initializing Ahats
    Ahat_hf{j} = zeros(r, r);
    Ahat_mf{j} = zeros(r, r);

    % Single fi Data
    X = Xhr(:, 1:m(j)*timesteps);
    Y = Xhrdot(:, 1:m(j)*timesteps)';

    % MFMC: Hifi Data 
    Xh = Xhr(:, 1:m1(j)*timesteps);
    Yh = Xhrdot(:, 1:m1(j)*timesteps)';

    % MFMC: Lofi Data
    Xl1 = Xlr(:, 1:m1(j)*timesteps);
    Yl1 = Xlrdot(:, 1:m1(j)*timesteps)';
    Xl2 = Xlr(:, 1:m2(j)*timesteps);
    Yl2 = Xlrdot(:, 1:m2(j)*timesteps)';
    
    for i = 1:r
        
        % single fidelity
        beta_single = XXT\((1./size(Y, 1))*X*Y(:, i));
        Ahat_hf{j}(i, :) = beta_single';

        % Multifidelity
        sigma1 = std(Xhr*Xhrdot(i, :)');
        sigma2 = std(Xlr*Xlrdot(i, :)');
        alpha = rhos(i)*sigma1./sigma2;
        XYlm1 = (1./size(Yl1, 1))*Xl1*Yl1(:, i);
        XYlm2 = (1./size(Yl2, 1))*Xl2*Yl2(:, i);
        XYhm1 = (1./size(Yh, 1))*Xh*Yh(:, i);
        XYmulti = XYhm1 + alpha*(XYlm2 - XYlm1);
        beta_multi = XXT\(XYmulti);
        Ahat_mf{j}(i, :) = beta_multi';
    end
    hf_errors(j) = norm(Ahat_hf{j} - Ahat_true, "fro")./norm(Ahat_true, "fro");
    mf_errors(j) = norm(Ahat_mf{j} - Ahat_true, "fro")./norm(Ahat_true, "fro");
end 
subplot(2, 4, r-1)
loglog(p, hf_errors, "o-")
hold on
loglog(p, mf_errors, "+-")
xlabel("Budget")
ylabel("Relative Error")
sgtitle("Budget versus Operator Error Relative to ")
title("Basis Size $r = " + r + "$")
legend("Single Fidelity", "Multifidelity")
xlim([min(p), max(p)])
end

