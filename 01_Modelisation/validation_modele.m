% validation_modele.m
% Simulation en boucle ouverte pour vérifier le comportement

clear; close all; clc;

% Charger les paramètres
run('parametres_systeme.m');
params = struct('M',M,'m',m,'L',L,'I',I,'b',b,'g',g);

% Conditions initiales (légèrement décalé de l'équilibre)
X0 = [0; 0.1; 0; 0];  % theta = 0.1 rad

% Temps de simulation
tspan = [0 5];

% Fonction d'entrée (force nulle)
F = 0;

% Simulation avec ode45
[t, X] = ode45(@(t,X) equations_mouvement(t, X, F, params), tspan, X0);

% Tracé
figure('Name','Boucle ouverte');
subplot(2,1,1);
plot(t, X(:,2)*180/pi); % theta en degrés
xlabel('Temps (s)');
ylabel('Angle (deg)');
title('Réponse en boucle ouverte à un écart initial');
grid on;

subplot(2,1,2);
plot(t, X(:,1));
xlabel('Temps (s)');
ylabel('Position chariot (m)');
grid on;

% Le système doit diverger (angle augmente)