% sim_LQG.m
% Simulation du pendule avec correcteur LQG (observateur par placement de pôles)

clear; close all; clc;

addpath('../../01_Modelisation');
run('../../01_Modelisation/parametres_systeme.m');
params = struct('M',M,'m',m,'L',L,'I',I,'b',b,'g',g);

load('gains_LQG.mat', 'K', 'L');

% Conditions initiales
X0_reel = [0; 0.1; 0; 0];
X0_est = [0; 0; 0; 0];

dt = 0.01;
t = 0:dt:10;
N = length(t);

X_reel = zeros(N,4);
X_est = zeros(N,4);
F = zeros(N,1);

X_reel(1,:) = X0_reel';
X_est(1,:) = X0_est';

C = [1 0 0 0; 0 1 0 0];
load('../../01_Modelisation/modele_lineaire.mat', 'A_num', 'B_num');

for i = 1:N-1
    F(i) = -K * X_est(i,:)';
    y = C * X_reel(i,:)' + 0.01 * randn(2,1);
    
    % Réel non linéaire
    dX_reel = equations_mouvement(t(i), X_reel(i,:)', F(i), params);
    X_reel(i+1,:) = X_reel(i,:) + dX_reel' * dt;
    
    % Estimation (observateur)
    dX_est = A_num * X_est(i,:)' + B_num * F(i) + L * (y - C * X_est(i,:)');
    X_est(i+1,:) = X_est(i,:) + dX_est' * dt;
end

% Tracés
figure('Name','LQG avec placement de pôles');
subplot(2,2,1);
plot(t, X_reel(:,2)*180/pi, 'b', t, X_est(:,2)*180/pi, 'r--');
xlabel('Temps (s)'); ylabel('Angle (deg)');
title('Angle réel vs estimé'); legend('Réel','Estimé'); grid on;
ylim([-10 30]);

subplot(2,2,2);
plot(t, X_reel(:,1), 'b', t, X_est(:,1), 'r--');
xlabel('Temps (s)'); ylabel('Position (m)');
title('Position chariot'); grid on;

subplot(2,2,3);
plot(t, X_reel(:,4)*180/pi, 'b', t, X_est(:,4)*180/pi, 'r--');
xlabel('Temps (s)'); ylabel('Vitesse angulaire (deg/s)');
title('Vitesse angulaire'); grid on;

subplot(2,2,4);
plot(t, F, 'g');
xlabel('Temps (s)'); ylabel('Force (N)');
title('Commande LQG'); grid on;