% linearisation.m
% Calcule les matrices A, B du système linéarisé autour de θ=0

clear; close all; clc;

% Charger les paramètres
run('parametres_systeme.m');

% Définition des variables symboliques
syms x theta x_dot theta_dot F real
syms M_sym m_sym L_sym I_sym b_sym g_sym real

% Équations non linéaires (d'après les équations de Lagrange)
% Matrice d'inertie
D = [M_sym + m_sym, m_sym * L_sym * cos(theta);
     m_sym * L_sym * cos(theta), I_sym + m_sym * L_sym^2];

% Vecteur des forces
C = [F - b_sym * x_dot + m_sym * L_sym * theta_dot^2 * sin(theta);
     m_sym * g_sym * L_sym * sin(theta)];

% Résolution
acc = D \ C;
x_ddot = acc(1);
theta_ddot = acc(2);

% Vecteur d'état
X = [x; theta; x_dot; theta_dot];
U = F;

% Dérivées
dX = [x_dot; theta_dot; x_ddot; theta_ddot];

% Calcul des Jacobiennes (linéarisation)
A_sym = jacobian(dX, X);
B_sym = jacobian(dX, U);

% Substitution numérique
% On crée une liste de substitutions
subs_list = [M_sym, m_sym, L_sym, I_sym, b_sym, g_sym, x, theta, x_dot, theta_dot, F];
vals_list = [M, m, L, I, b, g, 0, 0, 0, 0, 0];

A_num = double(subs(A_sym, subs_list, vals_list));
B_num = double(subs(B_sym, subs_list, vals_list));

% Affichage
disp('Matrice A :');
disp(A_num);
disp('Matrice B :');
disp(B_num);

% Pôles
poles = eig(A_num);
disp('Pôles du système linéarisé :');
disp(poles);

% Sauvegarde
save('modele_lineaire.mat', 'A_num', 'B_num');

% Vérification de la commandabilité
CO = ctrb(A_num, B_num);
rang = rank(CO);
fprintf('Rang de la matrice de commandabilité : %d (sur %d)\n', rang, size(A_num,1));
if rang == size(A_num,1)
    disp('✅ Système commandable');
else
    disp('⚠️ Système non commandable');
end

% Vérification de l'observabilité (pour LQG plus tard)
O = obsv(A_num, eye(4));
rang_obs = rank(O);
fprintf('Rang de la matrice d''observabilité : %d (sur %d)\n', rang_obs, size(A_num,1));