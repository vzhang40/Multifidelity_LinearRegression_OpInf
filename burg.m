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

t = linspace(0, T, K(1)+1);

xi = getRandIC(60, 1);
x0 = getIC(x, xi);
Xh = semiEuler(operators.A, operators.F,  x0, t);

i = 12;
Ur = U(:, 1:i);
x0r = Ur'*x0;
Oh = Ohat_hf{5};
Om = Ohat_mf{5};
esti_rstateHF = semiEuler(Oh(:, 1:i),Oh(:, i+1:end),  x0r, t);
esti_rstateMF = semiEuler(Om(:, 1:i),Om(:, i+1:end),  x0r, t);

sH = Ur*esti_rstateHF;
sM = Ur*esti_rstateMF;


%% ---- Fine Mesh Heatmap ----
figure
subplot(2,3,1)
colormap parula
imagesc(t(1:2000), x, Xh(:, 1:2000))
set(gca, 'XTickLabel', [], 'YTickLabel', [])
colormap hot
colorbar

title("Full State")

xlabel('$t$')

ylabel('$\eta$')
axis tight
axis equal

subplot(2,3,2)
imagesc(t(1:2000), x, sH(:, 1:2000))
set(gca, 'XTickLabel', [], 'YTickLabel', [])
colormap 
colorbar

title("Single Fidelity Operator Inference")

xlabel('$t$')

ylabel('$\eta$')
axis tight
axis equal

subplot(2,3,3)
imagesc(t(1:2000), x, sM(:, 1:2000))
set(gca, 'XTickLabel', [], 'YTickLabel', [])
colormap parula
colorbar

title("Multifidelity Operator Inference")

xlabel('$t$')

ylabel('$\eta$')
axis tight
axis equal

subplot(2,3,5)
imagesc(t(1:2000), x, Xh(:, 1:2000) - sH(:, 1:2000))
set(gca, 'XTickLabel', [], 'YTickLabel', [])

colorbar 
m = 0.02

caxis([-m m])

title("Single Fidelity Operator Inference Error")

xlabel('$t$')

ylabel('$\eta$')
axis tight
axis equal

subplot(2,3,6)
imagesc(t(1:2000), x, Xh(:, 1:2000) - sM(:, 1:2000))
set(gca, 'XTickLabel', [], 'YTickLabel', [])
colorbar
caxis([-m m])

title("Multifidelity Operator Inference Error")

xlabel('$t$')

ylabel('$\eta$')
axis tight
axis equal



%%    
    eH = 0.001*norm(Xh - sH, "fro")./norm(Xh, "fro");
    eM = 0.001*norm(Xh - sM, "fro")./norm(Xh, "fro");

%%
figure
 loglog(p, mean(eH, 2, 'omitnan'), "o-", "Color", rbg(1, :))
hold on
loglog(p, mean(eM, 2, 'omitnan'), "o-", "Color", rbg(2, :))
loglog(p,eH, "-", "Color", [rbg(1, :) 0.1], "LineWidth", 0.1)
loglog(p,eM, "-", "Color", [rbg(2, :) 0.1], "LineWidth", 0.1)
legend("Single Fidelity", "Multifidelity")

%%

m = max(abs(eH(:)));   % center around zero
caxis([-m m])

colormap(blueGrayOrangeBold(256))
colorbar

%%
function cmap = blueGrayOrangeBold(n)

% Bold diverging colormap: blue → gray → orange

if nargin < 1
    n = 256;
end

n1 = floor(n/2);
n2 = n - n1;

% Stronger endpoints
blue   = [0.10 0.30 0.85];
gray   = [0.92 0.92 0.92];
orange = [0.95 0.45 0.10];

% Negative side (blue → gray)
neg = [linspace(blue(1),gray(1),n1)', ...
       linspace(blue(2),gray(2),n1)', ...
       linspace(blue(3),gray(3),n1)'];

% Positive side (gray → orange)
pos = [linspace(gray(1),orange(1),n2)', ...
       linspace(gray(2),orange(2),n2)', ...
       linspace(gray(3),orange(3),n2)'];

cmap = [neg; pos];
end