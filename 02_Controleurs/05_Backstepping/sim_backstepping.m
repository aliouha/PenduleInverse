% sim_backstepping.m
% Simulation pendule inversé - Contrôle Backstepping
% Lancer depuis : 02_Controleurs/05_Backstepping/

clear; close all; clc;

% ── Chemins ────────────────────────────────────────────────────────────────
base = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(base);
addpath(fullfile(root, '01_Modelisation'));
addpath(fullfile(base, '05_Backstepping'));

% ── Paramètres physiques ───────────────────────────────────────────────────
run(fullfile(root, '01_Modelisation', 'parametres_systeme.m'));
params.M=M; params.m=m; params.L=L; params.I=I; params.b=b; params.g=g;

% ── Gains LQR (surface backstepping = surface LQR) ────────────────────────
load(fullfile(base, '02_LQR', 'gains_LQR.mat'), 'K');

% ── Paramètres Backstepping ────────────────────────────────────────────────
bs_params.K_lqr = K;   % surface s = K*X
bs_params.kb    = 5;   % gain convergence de s (s_dot = -kb*s)
                       % kb plus grand = convergence plus rapide
                       % trop grand = risque de saturation

% ── Condition initiale ─────────────────────────────────────────────────────
X0 = [0; 10*pi/180; 0; 0];

% ── Simulation ────────────────────────────────────────────────────────────
t_span  = [0, 10];
options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.01);

[t, X] = ode45(@(t,Xv) equations_mouvement(t, Xv, ...
    ctrl_backstepping(Xv, params, bs_params), params), ...
    t_span, X0, options);

% ── Recalcul commande et surface ───────────────────────────────────────────
F_hist = zeros(length(t),1);
s_hist = zeros(length(t),1);
for i = 1:length(t)
    F_hist(i) = ctrl_backstepping(X(i,:)', params, bs_params);
    s_hist(i) = bs_params.K_lqr * X(i,:)';
end

% ── Métriques ──────────────────────────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('   RESULTATS BACKSTEPPING\n');
fprintf('========================================\n');

idx1 = find(abs(X(:,2)) < 1*pi/180, 1, 'first');
if ~isempty(idx1)
    fprintf('Stabilisation theta < 1 deg   : %.2f s\n', t(idx1));
else
    fprintf('theta non stabilise a 1 deg\n');
end

idx_x = find(abs(X(:,1)) < 0.01, 1, 'first');
if ~isempty(idx_x)
    fprintf('Stabilisation x < 1 cm        : %.2f s\n', t(idx_x));
else
    fprintf('x non stabilise a 1 cm\n');
end

fprintf('Deplacement max chariot : %.1f cm\n', max(abs(X(:,1)))*100);
fprintf('Force max               : %.2f N\n',  max(abs(F_hist)));
fprintf('Chattering F fin        : %.5f N\n',  std(F_hist(end-50:end)));

if max(abs(X(:,2))) > 90*pi/180
    fprintf('ATTENTION : pendule a chute !\n');
else
    fprintf('Angle max : %.2f deg\n', max(abs(X(:,2)))*180/pi);
end
fprintf('========================================\n\n');

% ── Figures ────────────────────────────────────────────────────────────────
figure('Name','Backstepping','NumberTitle','off','Position',[50 50 1300 700]);

subplot(2,3,1)
plot(t, X(:,2)*180/pi, 'r', 'LineWidth',2);
yline(0,'k--'); yline(1,'b:','1 deg'); yline(-1,'b:');
xlabel('Temps (s)'); ylabel('theta (deg)');
title('Angle du pendule'); grid on;

subplot(2,3,2)
plot(t, X(:,1)*100, 'b', 'LineWidth',2);
yline(0,'k--');
xlabel('Temps (s)'); ylabel('x (cm)');
title('Position du chariot'); grid on;

subplot(2,3,3)
plot(t, F_hist, 'm', 'LineWidth',1.5);
yline(20,'r:'); yline(-20,'r:');
xlabel('Temps (s)'); ylabel('F (N)');
title('Commande (force)'); grid on;

subplot(2,3,4)
plot(t, s_hist, 'k', 'LineWidth',1.5);
yline(0,'r--','s=0');
xlabel('Temps (s)'); ylabel('s(t)');
title('Surface (doit converger vers 0)'); grid on;

subplot(2,3,5)
plot(t, X(:,4)*180/pi, 'r', 'LineWidth',1.5); hold on;
plot(t, X(:,3)*100,    'b', 'LineWidth',1.5);
legend('thetadot (deg/s)','xdot (cm/s)','Location','best');
xlabel('Temps (s)');
title('Vitesses'); grid on;

subplot(2,3,6)
angles_test = [5, 10, 15, 20];
colors_ci   = {'b','r','g','m'};
hold on;
for k = 1:length(angles_test)
    X0k = [0; angles_test(k)*pi/180; 0; 0];
    try
        [tk,Xk] = ode45(@(t,Xv) equations_mouvement(t,Xv,...
            ctrl_backstepping(Xv,params,bs_params),params),...
            [0,8], X0k, options);
        plot(tk, Xk(:,2)*180/pi, colors_ci{k}, 'LineWidth',1.5,...
            'DisplayName', sprintf('theta0=%d deg', angles_test(k)));
    catch
        fprintf('Echec theta0=%d deg\n', angles_test(k));
    end
end
yline(0,'k--');
xlabel('Temps (s)'); ylabel('theta (deg)');
legend('Location','best');
title('Robustesse - differentes CI'); grid on;

sgtitle('Controle par Backstepping','FontSize',14,'FontWeight','bold');