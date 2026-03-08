% calcul_LQR.m
% Calcul du gain optimal LQR pour le pendule inversé

clear; close all; clc;

% Ajouter le chemin
addpath('../../01_Modelisation');

% Charger paramètres
run('../../01_Modelisation/parametres_systeme.m');

% Charger le modèle linéarisé
load('../../01_Modelisation/modele_lineaire.mat', 'A_num', 'B_num');

% Définir les matrices de poids Q et R
% Q diagonale : poids sur [x, theta, x_dot, theta_dot]
% Plus le poids est grand, plus on pénalise cet état

% Version 1 : équilibrée (inspirée du PID réglé)
Q = diag([10, 200, 1, 10]);  % Poids fort sur theta (200), moyen sur x (10)
R = 0.1;                      % Poids sur l'effort de commande

% Calcul du gain LQR
[K, S, poles] = lqr(A_num, B_num, Q, R);

% Affichage
disp('Gain LQR K :');
disp(K);
disp('Pôles en boucle fermée :');
disp(poles);

% Sauvegarde
save('gains_LQR.mat', 'K', 'Q', 'R');

% Vérification de la stabilité
A_bf = A_num - B_num * K;
poles_bf = eig(A_bf);
if all(real(poles_bf) < 0)
    disp('✅ Système en boucle fermée stable');
else
    disp('⚠️ Instabilité !');
end