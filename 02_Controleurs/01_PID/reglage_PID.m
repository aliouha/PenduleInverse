% reglage_PID.m
% Gains PID corrigés pour le pendule inversé
% Basé sur l'analyse des pôles LQR : poles à -1, -4.5, -10.2

clear; close all; clc;

% ── Gains PID ─────────────────────────────────────────────────────────────
%
% Stratégie : deux boucles cascadées
%   Boucle interne (rapide) : stabilise theta
%   Boucle externe (lente)  : ramène x à 0 via une consigne theta
%
% F = Kp_theta*(0-theta) + Kd_theta*(0-theta_dot)   ← boucle theta
%   + Kp_x*(0-x) + Kd_x*(0-x_dot)                  ← boucle x
%
% Règle de réglage :
%   - Kp_theta doit compenser la gravité : m*g*L/delta ~ 3.8 rad/s²
%     → Kp_theta > (M+m)*delta/... → commencer à 100-200
%   - Kd_theta ~ 2*sqrt(Kp_theta) pour amortissement critique
%   - Kp_x ~ Kp_theta/10 (boucle lente)
%   - Kd_x ~ 2*sqrt(Kp_x)

Kp_theta = 200;    % gain proportionnel angle  (augmenté)
Kd_theta = 40;     % gain dérivé angle         (amorti)
Ki_theta = 5;      % gain intégral angle        (faible, anti-dérive)
Kp_x     = 10;     % gain proportionnel position
Kd_x     = 8;      % gain dérivé position       (NOUVEAU - amortit x_dot)
F_max    = 30;     % saturation [N]

save('gains_PID.mat', 'Kp_theta','Kd_theta','Ki_theta','Kp_x','Kd_x','F_max');

fprintf('Gains PID sauvegardes\n');
fprintf('  Kp_theta = %g\n', Kp_theta);
fprintf('  Kd_theta = %g\n', Kd_theta);
fprintf('  Ki_theta = %g\n', Ki_theta);
fprintf('  Kp_x     = %g\n', Kp_x);
fprintf('  Kd_x     = %g\n', Kd_x);
fprintf('  F_max    = %g N\n', F_max);