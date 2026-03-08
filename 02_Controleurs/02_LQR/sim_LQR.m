% sim_LQR.m
% Simulation du pendule avec correcteur LQR (retour d'état)

clear; close all; clc;

% Ajouter les chemins
addpath('../../01_Modelisation');

% Charger paramètres
run('../../01_Modelisation/parametres_systeme.m');
params = struct('M',M,'m',m,'L',L,'I',I,'b',b,'g',g);

% Charger le gain LQR
load('gains_LQR.mat', 'K');

% Conditions initiales (mêmes que pour le PID)
X0 = [0; 0.1; 0; 0];  % angle initial 0.1 rad (~5.7°)

% Temps
tspan = [0 10];

% Simulation avec ode45
[t, X] = ode45(@(t,X) pendule_LQR(t, X, params, K), tspan, X0);

% Calcul de la commande pour tracer
F = zeros(length(t),1);
for i = 1:length(t)
    F(i) = commande_LQR(X(i,:)', K);
end

% Tracé
figure('Name','LQR');
subplot(3,1,1);
plot(t, X(:,2)*180/pi, 'b', 'LineWidth', 1.5);
xlabel('Temps (s)'); ylabel('Angle (deg)');
title('Angle du pendule - LQR');
grid on;

subplot(3,1,2);
plot(t, X(:,1), 'r', 'LineWidth', 1.5);
xlabel('Temps (s)'); ylabel('Position (m)');
title('Position du chariot - LQR');
grid on;

subplot(3,1,3);
plot(t, F, 'g', 'LineWidth', 1.5);
xlabel('Temps (s)'); ylabel('Force (N)');
title('Commande - LQR');
grid on;

% Métriques de performance
[~, idx_max] = max(abs(X(:,2)));
depassement = (max(abs(X(:,2))) - abs(X(end,2))) * 180/pi;
idx_reponse = find(abs(X(:,2)*180/pi) < 2, 1, 'first');
if ~isempty(idx_reponse)
    temps_reponse = t(idx_reponse);
else
    temps_reponse = NaN;
end

fprintf('\n--- PERFORMANCES LQR ---\n');
fprintf('Dépassement : %.2f deg\n', depassement);
fprintf('Temps de réponse (seuil 2°) : %.2f s\n', temps_reponse);
fprintf('Erreur stationnaire : %.3f rad\n', abs(X(end,2)));

% Comparaison avec PID (si disponible)
try
    load('../01_PID/gains_PID.mat');
    fprintf('\n📊 Comparaison avec PID (valeurs finales) :\n');
    fprintf('   PID Kp_theta = %.1f, Kd_theta = %.1f\n', Kp_theta, Kd_theta);
    fprintf('   LQR gains : [%.1f, %.1f, %.1f, %.1f]\n', K(1), K(2), K(3), K(4));
catch
    disp('(Pas de données PID pour comparaison)');
end

% Fonction de commande LQR
function F = commande_LQR(X, K)
    % Retour d'état : u = -K * x
    F = -K * X;
end

% Fonction pour ode45
function dX = pendule_LQR(t, X, params, K)
    F = commande_LQR(X, K);
    dX = equations_mouvement(t, X, F, params);
end