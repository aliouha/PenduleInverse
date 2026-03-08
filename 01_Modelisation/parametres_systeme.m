% parametres_systeme.m
% Paramètres physiques du pendule inversé sur chariot

% Masses (kg)
M = 1.5;      % masse du chariot
m = 0.5;      % masse du pendule

% Longueurs (m)
L = 0.5;      % distance axe – centre de masse du pendule

% Inertie du pendule (kg*m^2) – si tige fine, I = m*L^2/3
I = m * L^2 / 3;

% Frottement (N/(m/s)) – coefficient de frottement visqueux du chariot
b = 0.1;

% Gravité (m/s^2)
g = 9.81;

% Paramètres pour les contrôleurs (seront ajustés plus tard)
% Q et R pour LQR
Q = diag([10, 10, 1, 1]);   % poids sur [x, theta, dx, dtheta]
R = 0.1;                     % poids sur la force

% Sauvegarde pour les autres scripts
save('parametres.mat', 'M', 'm', 'L', 'I', 'b', 'g', 'Q', 'R');

disp('Paramètres sauvegardés dans parametres.mat');