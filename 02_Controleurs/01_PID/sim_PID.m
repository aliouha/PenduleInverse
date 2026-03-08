% sim_PID.m
% Simulation pendule inversé - Contrôle PID
% Lancer depuis : 02_Controleurs/01_PID/

clear; close all; clc;

root = 'E:\Mes-CD\Matlab-Projects\PenduleInverse\PenduleInverse';
addpath(fullfile(root, '01_Modelisation'));
run(fullfile(root, '01_Modelisation', 'parametres_systeme.m'));
params = struct('M',M,'m',m,'L',L,'I',I,'b',b,'g',g);

load(fullfile(fileparts(mfilename('fullpath')), 'gains_PID.mat'));
gains = struct('Kp_theta',Kp_theta,'Kd_theta',Kd_theta,'Ki_theta',Ki_theta,...
               'Kp_x',Kp_x,'Kd_x',Kd_x,'F_max',F_max);

% ── Conditions initiales ──────────────────────────────────────────────────
X0 = [0; 0.1; 0; 0];   % 5.7°

% ── Simulation avec état intégral ─────────────────────────────────────────
% État étendu : [x, theta, x_dot, theta_dot, int_theta]
X0_ext = [X0; 0];   % intégrale theta = 0

tspan  = [0 10];
opts   = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.01);

[t, X_ext] = ode45(@(t,Xe) pendule_PID_ext(t, Xe, params, gains), ...
                   tspan, X0_ext, opts);

X     = X_ext(:,1:4);
F_cmd = zeros(length(t),1);
for i = 1:length(t)
    F_cmd(i) = commande_PID_ext(X_ext(i,:)', gains);
end

% ── Métriques ─────────────────────────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('   RESULTATS PID\n');
fprintf('========================================\n');

idx1 = find(abs(X(:,2)) < 1*pi/180, 1, 'first');
if ~isempty(idx1)
    fprintf('Stabilisation theta < 1 deg : %.2f s\n', t(idx1));
else
    fprintf('theta non stabilise\n');
end
fprintf('Deplacement max chariot : %.1f cm\n', max(abs(X(:,1)))*100);
fprintf('Force max               : %.2f N\n',  max(abs(F_cmd)));
if max(abs(X(:,2))) > 90*pi/180
    fprintf('ATTENTION : pendule a chute !\n');
else
    fprintf('Angle max : %.2f deg\n', max(abs(X(:,2)))*180/pi);
end
fprintf('========================================\n\n');

% ── Figures ───────────────────────────────────────────────────────────────
figure('Name','PID','NumberTitle','off','Position',[100 100 1000 700]);

subplot(2,2,1)
plot(t, X(:,2)*180/pi, 'b', 'LineWidth',2);
yline(0,'k--'); yline(1,'b:'); yline(-1,'b:');
xlabel('Temps (s)'); ylabel('theta (deg)');
title('Angle du pendule'); grid on;

subplot(2,2,2)
plot(t, X(:,1)*100, 'r', 'LineWidth',2);
yline(0,'k--');
xlabel('Temps (s)'); ylabel('x (cm)');
title('Position du chariot'); grid on;

subplot(2,2,3)
plot(t, F_cmd, 'g', 'LineWidth',1.5);
yline(gains.F_max,'r:'); yline(-gains.F_max,'r:');
xlabel('Temps (s)'); ylabel('F (N)');
title('Commande PID'); grid on;

subplot(2,2,4)
plot(t, X(:,4)*180/pi, 'r', 'LineWidth',1.5); hold on;
plot(t, X(:,3)*100,    'b', 'LineWidth',1.5);
legend('thetadot (deg/s)','xdot (cm/s)','Location','best');
xlabel('Temps (s)'); title('Vitesses'); grid on;

sgtitle('Controle PID - Pendule Inverse','FontSize',13,'FontWeight','bold');

% ── Fonctions locales ──────────────────────────────────────────────────────
function F = commande_PID_ext(Xe, gains)
    x         = Xe(1);
    theta     = Xe(2);
    x_dot     = Xe(3);
    theta_dot = Xe(4);
    int_theta = Xe(5);

    % Boucle theta (interne, rapide)
    F_theta = gains.Kp_theta * (0 - theta) ...
            + gains.Kd_theta * (0 - theta_dot) ...
            + gains.Ki_theta * int_theta;

    % Boucle x (externe, lente)
    F_x = gains.Kp_x * (0 - x) ...
        + gains.Kd_x * (0 - x_dot);

    F = F_theta + F_x;
    F = max(min(F, gains.F_max), -gains.F_max);
end

function dXe = pendule_PID_ext(t, Xe, params, gains)
    F  = commande_PID_ext(Xe, gains);
    dX = equations_mouvement(t, Xe(1:4), F, params);
    % Intégrale de theta
    d_int_theta = Xe(2);   % d(int_theta)/dt = theta
    dXe = [dX; d_int_theta];
end