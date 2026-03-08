% calcul_LQG.m
% Synthèse du correcteur LQG avec observateur par placement de pôles

clear; close all; clc;

% Ajouter les chemins
addpath('../../01_Modelisation');

% Charger paramètres
run('../../01_Modelisation/parametres_systeme.m');

% Charger le modèle linéarisé
load('../../01_Modelisation/modele_lineaire.mat', 'A_num', 'B_num');

% ==================== LQR ====================
Q_lqr = diag([10, 200, 1, 10]);
R_lqr = 0.1;
[K, ~, ~] = lqr(A_num, B_num, Q_lqr, R_lqr);
disp('Gain LQR K :');
disp(K);

% ==================== Observateur par placement de pôles ====================
% On mesure la position du chariot et l'angle
C = [1 0 0 0; 
     0 1 0 0];

% Choix des pôles de l'observateur
% On les choisit 3 à 5 fois plus rapides que les pôles du LQR
poles_lqr = eig(A_num - B_num*K);
disp('Pôles du LQR :');
disp(poles_lqr);

% Pôles de l'observateur (tous réels négatifs, plus rapides)
poles_obs = [-8, -9, -10, -11];  % À ajuster si nécessaire

% Calcul du gain par placement de pôles
L = place(A_num', C', poles_obs)';
disp('Gain de l''observateur L :');
disp(L);

% Vérification de la stabilité
if all(real(eig(A_num - L*C)) < 0)
    disp('✅ Observateur stable');
else
    disp('❌ Observateur instable');
end

% Sauvegarde
save('gains_LQG.mat', 'K', 'L');