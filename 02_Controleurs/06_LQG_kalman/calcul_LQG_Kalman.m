% calcul_LQG_Kalman.m
% Synthèse LQG avec filtre de Kalman (observateur optimal stochastique)
% Différence avec placement de pôles :
%   - Luenberger : L choisi arbitrairement → pas optimal
%   - Kalman     : L calculé pour minimiser E[||X - X_hat||²] → optimal
%
% Modèle stochastique :
%   x_dot = A*x + B*u + G*w    (w = bruit de process,  cov = Qk)
%   y     = C*x + v            (v = bruit de mesure,   cov = Rk)

clear; close all; clc;

base = fileparts(mfilename('fullpath'));   % → 06_LQG_kalman
root = 'E:\Mes-CD\Matlab-Projects\PenduleInverse\PenduleInverse';

addpath(fullfile(root, '01_Modelisation'));
run(fullfile(root, '01_Modelisation', 'parametres_systeme.m'));
load(fullfile(root, '01_Modelisation', 'modele_lineaire.mat'), 'A_num', 'B_num');

% ── 1. Correcteur LQR (identique à avant) ─────────────────────────────────
Q_lqr = diag([10, 200, 1, 10]);
R_lqr = 0.1;
[K, ~, ~] = lqr(A_num, B_num, Q_lqr, R_lqr);
fprintf('Gain LQR K :\n'); disp(K);

poles_lqr = eig(A_num - B_num*K);
fprintf('Poles LQR en BF :\n'); disp(poles_lqr);

% ── 2. Matrice d'observation (mesures disponibles) ────────────────────────
% On mesure x (position) et theta (angle) — mêmes capteurs que ton LQG
C = [1 0 0 0;   % mesure x
     0 1 0 0];  % mesure theta

% Vérification de l'observabilité
if rank(obsv(A_num, C)) == size(A_num,1)
    fprintf('Systeme observable avec C\n');
else
    fprintf('ATTENTION : systeme non observable !\n');
end

% ── 3. Paramètres du bruit (réglage Kalman) ───────────────────────────────
%
% Qk = covariance du bruit de process (perturbations sur le système)
%      Grande Qk → on fait moins confiance au modèle → Kalman plus réactif
%      Petite Qk → on fait confiance au modèle → Kalman plus lisse
%
% Rk = covariance du bruit de mesure (bruit capteurs)
%      Grande Rk → capteurs peu fiables → Kalman corrige peu
%      Petite Rk → capteurs précis → Kalman suit les mesures
%
% Règle pratique :
%   - Rk se lit directement sur les capteurs réels (variance du bruit)
%   - Qk se règle pour obtenir le compromis vitesse/lissage voulu
%   - Ratio Qk/Rk contrôle l'agressivité du filtre

% Bruit de process : perturbations sur chaque état
% [x, theta, x_dot, theta_dot]
Qk = diag([0.01,   ...  % incertitude sur x (faible)
           0.01,   ...  % incertitude sur theta (faible)
           0.1,    ...  % incertitude sur x_dot (plus grande)
           0.1]);       % incertitude sur theta_dot (plus grande)

% Bruit de mesure : bruit des capteurs (position et angle)
Rk = diag([0.001,  ...  % capteur position x    (précis)
           0.001]);     % capteur angle theta   (précis)

% ── 4. Calcul du gain de Kalman ───────────────────────────────────────────
%
% Le filtre de Kalman résout le problème dual du LQR :
%   min E[||X - X_hat||²]
%
% Mathématiquement : on résout l'équation de Riccati pour P (cov de l'erreur)
%   A*P + P*A' - P*C'*Rk^(-1)*C*P + Qk = 0
%   L = P * C' * Rk^(-1)
%
% En MATLAB, lqr() sur le système dual (A', C', Qk, Rk) donne L'
[~, P, ~] = lqr(A_num', C', Qk, Rk);
L_kalman = P * C' / Rk;

fprintf('\nGain de Kalman L :\n'); disp(L_kalman);

% ── 5. Comparaison avec le Luenberger ─────────────────────────────────────
poles_obs_luenberger = [-8, -9, -10, -11];
L_luenberger = place(A_num', C', poles_obs_luenberger)';

fprintf('Poles observateur Luenberger : [-8, -9, -10, -11]\n');
fprintf('Poles observateur Kalman     :\n');
disp(eig(A_num - L_kalman*C));

fprintf('Comparaison normes gains :\n');
fprintf('  ||L_luenberger|| = %.4f\n', norm(L_luenberger));
fprintf('  ||L_kalman||     = %.4f\n', norm(L_kalman));

% ── 6. Vérification stabilité observateur Kalman ──────────────────────────
poles_kalman = eig(A_num - L_kalman*C);
if all(real(poles_kalman) < 0)
    fprintf('Observateur Kalman stable\n');
else
    fprintf('ATTENTION : Observateur Kalman instable !\n');
end

% ── 7. Sauvegarde ─────────────────────────────────────────────────────────
save('gains_LQG_Kalman.mat', 'K', 'L_kalman', 'L_luenberger', ...
     'Qk', 'Rk', 'C', 'P');
fprintf('\nGains sauvegardes dans gains_LQG_Kalman.mat\n');