% comparaison_LQG.m
clear; close all; clc;

root = 'E:\Mes-CD\Matlab-Projects\PenduleInverse\PenduleInverse';
addpath(fullfile(root, '01_Modelisation'));
run(fullfile(root, '01_Modelisation', 'parametres_systeme.m'));
params = struct('M',M,'m',m,'L',L,'I',I,'b',b,'g',g);
load(fullfile(root, '01_Modelisation', 'modele_lineaire.mat'), 'A_num', 'B_num');
load(fullfile(fileparts(mfilename('fullpath')), 'gains_LQG_Kalman.mat'), ...
     'K', 'L_kalman', 'L_luenberger');

C  = [1 0 0 0; 0 1 0 0];
Rk = diag([0.001, 0.001]);
sigma_x      = sqrt(Rk(1,1));
sigma_theta  = sqrt(Rk(2,2));
sigma_process = 0;

dt = 0.005; T = 10; t = (0:dt:T)'; N = length(t);
X0_reel = [0; 0.1; 0; 0]; 

% Même bruit pour comparaison équitable
rng(42);
bruit_m = randn(N, 2);
bruit_p = randn(N, 1);

% ── Simulation Luenberger ─────────────────────────────────────────────────
Xr_lub = zeros(N,4); Xe_lub = zeros(N,4); F_lub = zeros(N,1);
Xr_lub(1,:) = X0_reel'; Xe_lub(1,:) = zeros(1,4);

for i = 1:N-1
    F_lub(i) = -K * Xe_lub(i,:)';
    y = C*Xr_lub(i,:)' + [sigma_x*bruit_m(i,1); sigma_theta*bruit_m(i,2)];
    dXr = equations_mouvement(t(i), Xr_lub(i,:)', F_lub(i)+sigma_process*bruit_p(i), params);
    Xr_lub(i+1,:) = Xr_lub(i,:) + dXr'*dt;
    inno = y - C*Xe_lub(i,:)';
    dXe = A_num*Xe_lub(i,:)' + B_num*F_lub(i) + L_luenberger*inno;
    Xe_lub(i+1,:) = Xe_lub(i,:) + dXe'*dt;
end

% ── Simulation Kalman ─────────────────────────────────────────────────────
Xr_kal = zeros(N,4); Xe_kal = zeros(N,4); F_kal = zeros(N,1);
Xr_kal(1,:) = X0_reel'; Xe_kal(1,:) = zeros(1,4);

for i = 1:N-1
    F_kal(i) = -K * Xe_kal(i,:)';
    y = C*Xr_kal(i,:)' + [sigma_x*bruit_m(i,1); sigma_theta*bruit_m(i,2)];
    dXr = equations_mouvement(t(i), Xr_lub(i,:)', F_lub(i), params);
    Xr_kal(i+1,:) = Xr_kal(i,:) + dXr'*dt;
    inno = y - C*Xe_kal(i,:)';
    dXe = A_num*Xe_kal(i,:)' + B_num*F_kal(i) + L_kalman*inno;
    Xe_kal(i+1,:) = Xe_kal(i,:) + dXe'*dt;
end

% ── Métriques ─────────────────────────────────────────────────────────────
E_lub = Xr_lub - Xe_lub;
E_kal = Xr_kal - Xe_kal;

fprintf('\n==================================================\n');
fprintf('              LUENBERGER    KALMAN\n');
fprintf('==================================================\n');
fprintf('Err theta RMS : %.4f deg    %.4f deg\n', ...
    rms(E_lub(:,2))*180/pi, rms(E_kal(:,2))*180/pi);
fprintf('Err x RMS     : %.4f cm     %.4f cm\n', ...
    rms(E_lub(:,1))*100, rms(E_kal(:,1))*100);
fprintf('Force max     : %.2f N       %.2f N\n', ...
    max(abs(F_lub)), max(abs(F_kal)));
fprintf('==================================================\n\n');

% ── Figures ───────────────────────────────────────────────────────────────
figure('Name','Comparaison LQG','Position',[50 50 1300 600]);

subplot(1,3,1)
plot(t, Xr_lub(:,2)*180/pi, 'b',  'LineWidth',2); hold on;
plot(t, Xr_kal(:,2)*180/pi, 'r--','LineWidth',2);
yline(0,'k--'); xlabel('Temps (s)'); ylabel('theta (deg)');
legend('Luenberger','Kalman','Location','best');
title('Angle pendule (etat reel)'); grid on;

subplot(1,3,2)
plot(t, Xr_lub(:,1)*100, 'b',  'LineWidth',2); hold on;
plot(t, Xr_kal(:,1)*100, 'r--','LineWidth',2);
yline(0,'k--'); xlabel('Temps (s)'); ylabel('x (cm)');
legend('Luenberger','Kalman');
title('Position chariot (etat reel)'); grid on;

subplot(1,3,3)
plot(t, E_lub(:,2)*180/pi, 'b',  'LineWidth',1.5); hold on;
plot(t, E_kal(:,2)*180/pi, 'r--','LineWidth',1.5);
yline(0,'k--'); xlabel('Temps (s)'); ylabel('erreur theta (deg)');
legend('Luenberger','Kalman');
title('Erreur estimation theta'); grid on;

sgtitle('Comparaison LQG : Luenberger vs Kalman', ...
        'FontSize',13,'FontWeight','bold');