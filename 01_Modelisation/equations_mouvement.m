function dX = equations_mouvement(t, X, F, params)
% EQUATIONS_MOUVEMENT  Modèle non linéaire du pendule inversé sur chariot
%   t : temps (inutilisé mais nécessaire pour ode45)
%   X : vecteur d'état [x; theta; x_dot; theta_dot]
%   F : force appliquée au chariot (entrée)
%   params : structure contenant M, m, L, I, b, g
%   dX : dérivées [x_dot; theta_dot; x_ddot; theta_ddot]

% Extraction des états
x = X(1);
theta = X(2);
x_dot = X(3);
theta_dot = X(4);

% Paramètres
M = params.M;
m = params.m;
L = params.L;
I = params.I;
b = params.b;
g = params.g;

% Calcul des accélérations (d'après les équations de Lagrange)
% On utilise les formules classiques (voir par exemple : 
% https://underactuated.mit.edu/pend.html)
% Pour simplifier, on utilise la forme matricielle :
% [M+m, m*L*cos(theta); m*L*cos(theta), I+m*L^2] * [x_ddot; theta_ddot] 
%   = [F - b*x_dot + m*L*theta_dot^2*sin(theta); 
%      m*g*L*sin(theta)]

% Matrice d'inertie
D = [M+m, m*L*cos(theta);
     m*L*cos(theta), I + m*L^2];

% Vecteur des forces
C = [F - b*x_dot + m*L*theta_dot^2*sin(theta);
     m*g*L*sin(theta)];

% Résolution pour les accélérations
acc = D \ C;

x_ddot = acc(1);
theta_ddot = acc(2);

% Dérivées
dX = [x_dot; theta_dot; x_ddot; theta_ddot];
end