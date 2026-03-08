% sim_LQG_Kalman.m
% Simulation LQG avec filtre de Kalman
% Compare : système réel vs état estimé par Kalman
% Lancer depuis : 02_Controleurs/03_LQG/

clear; close all; clc;

root = 'E:\Mes-CD\Matlab-Projects\PenduleInverse\PenduleInverse';
addpath(fullfile(root, '01_Modelisation'));
run(fullfile(root, '01_Modelisation', 'parametres_systeme.m'));
params = struct('M',M,'m',m,'L',L,'I',I,'b',b,'g',g);
load(fullfile(root, '01_Modelisation', 'modele_lineaire.mat'), 'A_num', 'B_num');
load(fullfile(fileparts(mfilename('fullpath')), 'gains_LQG_Kalman.mat'), 'K', 'L_kalman');
Rk = diag([0.001, 0.001]);
C  = [1 0 0 0; 0 1 0 0];

% ── Paramètres simulation ─────────────────────────────────────────────────
dt  = 0.005;          % pas de temps (plus petit = plus précis)
T   = 10;             % durée totale
t   = (0:dt:T)';
N   = length(t);

% ── Conditions initiales ──────────────────────────────────────────────────
X0_reel = [0; 0.1; 0; 0];    % état réel initial (0.1 rad ≈ 5.7°)
X0_est  = [0; 0;   0; 0];    % estimé initial (observateur part de 0)

% ── Niveaux de bruit ──────────────────────────────────────────────────────
% Bruit de mesure (écart-type = sqrt(Rk))
sigma_x     = sqrt(Rk(1,1));   % bruit sur x
sigma_theta = sqrt(Rk(2,2));   % bruit sur theta

% Bruit de process (perturbation externe sur le système)
sigma_process = 0.1;   % [N] perturbation force aléatoire

% ── Initialisation ────────────────────────────────────────────────────────
X_reel = zeros(N, 4);
X_est  = zeros(N, 4);
F_hist = zeros(N, 1);
Y_hist = zeros(N, 2);   % mesures bruitées
E_hist = zeros(N, 4);   % erreur d'estimation

X_reel(1,:) = X0_reel';
X_est(1,:)  = X0_est';

% ── Boucle de simulation ──────────────────────────────────────────────────
for i = 1:N-1

    % 1. Commande LQR sur état ESTIMÉ (pas le réel → c'est le LQG)
    F_hist(i) = -K * X_est(i,:)';

    % 2. Mesures bruitées (capteurs imparfaits)
    bruit_mesure = [sigma_x * randn;
                    sigma_theta * randn];
    Y_hist(i,:) = (C * X_reel(i,:)' + bruit_mesure)';

    % 3. Évolution du système RÉEL non linéaire + bruit de process
    bruit_process = sigma_process * randn;   % perturbation force
    dX_reel = equations_mouvement(t(i), X_reel(i,:)', ...
                                  F_hist(i) + bruit_process, params);
    X_reel(i+1,:) = X_reel(i,:) + dX_reel' * dt;

    % 4. Filtre de Kalman (modèle LINÉAIRE + correction mesures)
    %    x̂_dot = A*x̂ + B*u + L*(y - C*x̂)
    %    Innovation = y - C*x̂  (écart mesure / prédiction)
    innovation = Y_hist(i,:)' - C * X_est(i,:)';
    dX_est = A_num * X_est(i,:)' + B_num * F_hist(i) + L_kalman * innovation;
    X_est(i+1,:)  = X_est(i,:) + dX_est' * dt;

    % 5. Erreur d'estimation
    E_hist(i,:) = X_reel(i,:) - X_est(i,:);
end

% Dernière mesure
Y_hist(N,:) = (C * X_reel(N,:)')';
E_hist(N,:) = X_reel(N,:) - X_est(N,:);

% ── Métriques ─────────────────────────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('   RESULTATS LQG - FILTRE DE KALMAN\n');
fprintf('========================================\n');

idx1 = find(abs(X_reel(:,2)) < 1*pi/180, 1, 'first');
if ~isempty(idx1)
    fprintf('Stabilisation theta < 1 deg : %.2f s\n', t(idx1));
end

fprintf('Deplacement max chariot     : %.1f cm\n', max(abs(X_reel(:,1)))*100);
fprintf('Force max                   : %.2f N\n',  max(abs(F_hist)));
fprintf('Erreur estimation theta RMS : %.4f deg\n', ...
    rms(E_hist(:,2))*180/pi);
fprintf('Erreur estimation x RMS     : %.4f cm\n', ...
    rms(E_hist(:,1))*100);
fprintf('========================================\n\n');

% ── Figures ───────────────────────────────────────────────────────────────
figure('Name','LQG Kalman','NumberTitle','off','Position',[50 50 1300 800]);

% Angle theta
subplot(3,3,1)
plot(t, X_reel(:,2)*180/pi, 'b',  'LineWidth',2); hold on;
plot(t, X_est(:,2)*180/pi,  'r--','LineWidth',1.5);
yline(0,'k--');
xlabel('Temps (s)'); ylabel('theta (deg)');
title('Angle pendule : reel vs estime');
legend('Reel','Kalman','Location','best'); grid on;

% Position x
subplot(3,3,2)
plot(t, X_reel(:,1)*100, 'b',  'LineWidth',2); hold on;
plot(t, X_est(:,1)*100,  'r--','LineWidth',1.5);
yline(0,'k--');
xlabel('Temps (s)'); ylabel('x (cm)');
title('Position chariot : reel vs estime');
legend('Reel','Kalman','Location','best'); grid on;

% Vitesse angulaire
subplot(3,3,3)
plot(t, X_reel(:,4)*180/pi, 'b',  'LineWidth',2); hold on;
plot(t, X_est(:,4)*180/pi,  'r--','LineWidth',1.5);
xlabel('Temps (s)'); ylabel('thetadot (deg/s)');
title('Vitesse angulaire : reel vs estime');
legend('Reel','Kalman','Location','best'); grid on;

% Mesures bruitées vs réel
subplot(3,3,4)
plot(t, X_reel(:,2)*180/pi,  'b', 'LineWidth',2); hold on;
plot(t, Y_hist(:,2)*180/pi, 'g.', 'MarkerSize',3);
plot(t, X_est(:,2)*180/pi,  'r--','LineWidth',1.5);
xlabel('Temps (s)'); ylabel('theta (deg)');
title('Mesures bruitees vs Kalman');
legend('Reel','Mesure bruitee','Kalman','Location','best'); grid on;

% Erreur d'estimation theta
subplot(3,3,5)
plot(t, E_hist(:,2)*180/pi, 'r', 'LineWidth',1.5);
yline(0,'k--');
xlabel('Temps (s)'); ylabel('erreur theta (deg)');
title('Erreur estimation theta');
grid on;

% Erreur d'estimation x
subplot(3,3,6)
plot(t, E_hist(:,1)*100, 'b', 'LineWidth',1.5);
yline(0,'k--');
xlabel('Temps (s)'); ylabel('erreur x (cm)');
title('Erreur estimation x');
grid on;

% Commande
subplot(3,3,7)
plot(t, F_hist, 'm', 'LineWidth',1.5);
xlabel('Temps (s)'); ylabel('F (N)');
title('Commande LQG'); grid on;

% Innovation (signal clé du Kalman)
subplot(3,3,8)
innovation_theta = Y_hist(:,2) - X_est(:,2);   % y_theta - C_theta*x_hat
plot(t, innovation_theta*180/pi, 'k', 'LineWidth',1.2);
yline(0,'r--');
xlabel('Temps (s)'); ylabel('innovation (deg)');
title('Innovation Kalman sur theta (y - C*xhat)');
grid on;

% Vitesse chariot estimée (non mesurée !)
subplot(3,3,9)
plot(t, X_reel(:,3)*100, 'b',  'LineWidth',2); hold on;
plot(t, X_est(:,3)*100,  'r--','LineWidth',1.5);
xlabel('Temps (s)'); ylabel('xdot (cm/s)');
title('Vitesse chariot (NON mesuree, estimee par Kalman)');
legend('Reelle','Estimee','Location','best'); grid on;

sgtitle('LQG avec Filtre de Kalman','FontSize',14,'FontWeight','bold');