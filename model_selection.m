%% Problem Set-Up
clear; close all; clc

%% 1D Periodic Heat Equation
addpath("data_files")
load("BurgerEqData.mat")

% Getting POD data with hifi data
[U, ~, ~] = svd(datah, "econ");

% Budgets
p = [10, 50, 100, 1000, 5000];

% Costs
w = [1; 0.01]; %16^2/(128^2)]; % about 0.0156

%% Change r for debugging
r = 12;

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

%% Model Statistics
G11 = cell(1, r);
G12 = cell(1, r);
G22 = cell(1, r);
alpha = zeros(r, 1);
for i = 1:r
    ghigh = Xhr.*Xhrdot(i, :);
    glow = Xlr.*Xlrdot(i, :);

    ghighvar = ghigh - mean(ghigh);
    glowvar = glow - mean(glow);
    
    Gamhh = (ghighvar * ghighvar') / (size(glowvar,2) - 1);
    Gamhl = (ghighvar * glowvar') / (size(glowvar,2) - 1);
    Gamll = (glowvar * glowvar') / (size(glowvar,2) - 1);
    alpha(i) = trace(Gamhl)/trace(Gamll);

    G11{i} = Gamhh;
    G12{i} = Gamhl;
    G22{i} = Gamll;
end

%% MFMC Model Statistics
rhos = zeros(r, 3);
for i = 1:r
    temp = corrcoef(Xhrdot(i, :)', Xlrdot(i, :)');
    rhos(i, 1:end-1) = temp(1, :);
end
sigma1 = std(Xhrdot, 1, 2);
sigma2 = std(Xlrdot, 1, 2);

%% Find "optimal" sample size, this is MFMC method
m1 = zeros(r, length(p));
m2 = zeros(r, length(p));
for i = 1:r
    rs = sqrt(w(1)*(rhos(i, 1:end-1).^2 - rhos(i, 2:end).^2) ./ (w'*(1 - rhos(i, 2).^2)))';
    m1(i, :) = floor(p./(w'*rs));
    m2(i, :) = floor(m1(i, :).*rs(2));
end
m = floor(p./w(1));


%% Trace of covariance of MC estimated vector
inds = 1:r;
% Previous numerical experiment: arbitrary m1 and m2
% m1 = repmat(0.8*p, [r, 1]);
% m2 = (p - m1)./w(2);

TrC = zeros(r, length(p));
TrC_MF = zeros(r, length(p));
for i = 1:length(p)
    for j = 1:r
        TrC(j, i) = (1./m(i)).*trace(G11{j});
        TrC_MF(j, i) = (1./m1(j, i)).*trace(G11{j}) + (1./m1(j, i) - 1./m2(j, i)).*(alpha(j)^2*trace(G22{j}) - 2*alpha(j)*trace(G12{j}));
    end
end
mask = all(0 < TrC - TrC_MF, 2);
disp("Rank sizes " + num2str(inds(~mask)) + " exhibit higher MF trace of covariance")

%% MFMC Model Selection

% testing cost versus correlation
RHS = (rhos(:, 1:end-2).^2 - rhos(:, 2:end-1).^2) ./ (rhos(:, 2:end-1).^2 - rhos(:, 3:end).^2);
mask1 = w(1)./w(2:end) > RHS;

% testing normalized with cost variance p = 1, w1 = 1
v_MC = sigma1.^2;
v_MFMC = ( sigma1.^2 ) .* (sqrt(w(1).*(rhos(:, 1).^2 - rhos(:, 2).^2)) + sqrt(w(2).*(rhos(:, 2).^2 - rhos(:, 3).^2)) ).^2;
mask2 = 0 < (v_MC - v_MFMC);

if all(mask1)
    disp("Condition (3.12) is satisfied")
else 
    disp("Rank sizes " + num2str(inds(~mask1)) + " violate Condition (3.12)")
end

if all(mask2)
    disp("MFMC works for all rank sizes")
else 
    disp("Rank sizes " + num2str(inds(~mask2)) + " exhibit higher MFMC variance")
end

