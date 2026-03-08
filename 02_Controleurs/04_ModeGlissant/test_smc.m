% test_smc.m
% Script de diagnostic complet pour le Mode Glissant
% Lancer depuis : 02_Controleurs/04_ModeGlissant/

clear; close all; clc;
fprintf('========================================\n');
fprintf('   DIAGNOSTIC - Mode Glissant (SMC)\n');
fprintf('========================================\n\n');

% ── 1. Chemins ────────────────────────────────────────────────────────────
base = fileparts(fileparts(mfilename('fullpath')));  % → 02_Controleurs
root = fileparts(base);                               % → PenduleInverse

addpath(fullfile(root, '01_Modelisation'));
addpath(fullfile(base, '04_ModeGlissant'));

fprintf('[1] Chemins ajoutés :\n');
fprintf('    Modelisation : %s\n', fullfile(root,'01_Modelisation'));
fprintf('    SMC          : %s\n\n', fullfile(base,'04_ModeGlissant'));

% ── 2. Paramètres physiques ───────────────────────────────────────────────
fprintf('[2] Chargement parametres_systeme.m... ');
try
    run(fullfile(root, '01_Modelisation', 'parametres_systeme.m'));
    params.M=M; params.m=m; params.L=L; params.I=I; params.b=b; params.g=g;
    fprintf('OK\n');
    fprintf('    M=%.2f kg, m=%.2f kg, L=%.2f m, g=%.2f\n\n', M,m,L,g);
catch e
    fprintf('ERREUR : %s\n', e.message); return;
end

% ── 3. Chargement gains LQR ──────────────────────────────────────────────
fprintf('[3] Chargement gains LQR... ');
try
    load(fullfile(base, '02_LQR', 'gains_LQR.mat'), 'K');
    fprintf('OK  →  K = [%.2f, %.2f, %.2f, %.2f]\n\n', K(1),K(2),K(3),K(4));
catch e
    fprintf('ERREUR : %s\n', e.message); return;
end

% ── 4. Paramètres SMC ────────────────────────────────────────────────────
smc_params.K_lqr = K;    % surface s = K*X
smc_params.k     = 8;    % gain robustesse
smc_params.phi   = 0.3;  % couche limite anti-chattering

% ── 5. Test contrôleur seul ───────────────────────────────────────────────
fprintf('[4] Test ctrl_mode_glissant sur X0... ');
X_test = [0; 10*pi/180; 0; 0];
try
    F_test = ctrl_mode_glissant(X_test, params, smc_params);
    s_test = smc_params.K_lqr * X_test;
    fprintf('OK  →  F = %.4f N,  s = %.4f\n\n', F_test, s_test);
catch e
    fprintf('ERREUR : %s\n', e.message); return;
end

% ── 6. Test equations_mouvement ──────────────────────────────────────────
fprintf('[5] Test equations_mouvement... ');
try
    dX = equations_mouvement(0, X_test, F_test, params);
    fprintf('OK  →  dX = [%.4f, %.4f, %.4f, %.4f]\n\n', ...
        dX(1),dX(2),dX(3),dX(4));
catch e
    fprintf('ERREUR : %s\n', e.message); return;
end

% ── 7. Simulation ODE45 ──────────────────────────────────────────────────
fprintf('[6] Simulation ODE45 (10s)... ');
t_span  = [0, 10];
options = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.01);

try
    [t, X] = ode45(@(t,Xv) equations_mouvement(t, Xv, ...
        ctrl_mode_glissant(Xv, params, smc_params), params), ...
        t_span, X_test, options);
    fprintf('OK  →  %d points calculés\n\n', length(t));
catch e
    fprintf('ERREUR ODE45 : %s\n', e.message); return;
end

% ── 8. Recalcul commande + surface ───────────────────────────────────────
F_hist = zeros(length(t),1);
s_hist = zeros(length(t),1);

for i = 1:length(t)
    F_hist(i) = ctrl_mode_glissant(X(i,:)', params, smc_params);
    s_hist(i) = smc_params.K_lqr * X(i,:)';
end

% ── 9. Métriques ─────────────────────────────────────────────────────────
fprintf('========================================\n');
fprintf('   MÉTRIQUES DE PERFORMANCE\n');
fprintf('========================================\n');

idx1 = find(abs(X(:,2)) < 1*pi/180, 1, 'first');
if ~isempty(idx1)
    fprintf('✅ Stabilisation θ < 1°   : %.2f s\n', t(idx1));
else
    fprintf('⚠️  θ non stabilisé à 1°\n');
end

idx05 = find(abs(X(:,2)) < 0.5*pi/180, 1, 'first');
if ~isempty(idx05)
    fprintf('✅ Stabilisation θ < 0.5° : %.2f s\n', t(idx05));
end

fprintf('📏 Déplacement max chariot  : %.4f m  (%.1f cm)\n', ...
    max(abs(X(:,1))), max(abs(X(:,1)))*100);
fprintf('⚡ Force max appliquée      : %.2f N\n', max(abs(F_hist)));
fprintf('🔄 Chattering (std F fin)   : %.5f N\n', std(F_hist(end-50:end)));
fprintf('📐 Angle max atteint        : %.2f°\n', max(abs(X(:,2)))*180/pi);

if max(abs(X(:,2))) > 90*pi/180
    fprintf('❌ Pendule a chuté (θ > 90°)\n');
elseif max(abs(X(:,2))) > 30*pi/180
    fprintf('⚠️  Angle max important : %.1f°\n', max(abs(X(:,2)))*180/pi);
else
    fprintf('✅ Angle bien contenu\n');
end
fprintf('========================================\n\n');

% ── 10. Figures ───────────────────────────────────────────────────────────
figure('Name','SMC - Diagnostic','NumberTitle','off','Position',[50 50 1300 700]);

subplot(2,3,1)
plot(t, X(:,2)*180/pi, 'r', 'LineWidth',2);
yline(0,'k--','Réf'); yline(1,'b:','1°'); yline(-1,'b:');
xlabel('Temps (s)'); ylabel('\theta (°)');
title('Angle du pendule'); grid on;

subplot(2,3,2)
plot(t, X(:,1)*100, 'b', 'LineWidth',2);
yline(0,'k--');
xlabel('Temps (s)'); ylabel('x (cm)');
title('Position du chariot'); grid on;

subplot(2,3,3)
plot(t, F_hist, 'm', 'LineWidth',1.5);
yline( 20,'r:','+F_{max}'); yline(-20,'r:','-F_{max}');
xlabel('Temps (s)'); ylabel('F (N)');
title('Commande (force)'); grid on;

subplot(2,3,4)
plot(t, s_hist, 'k', 'LineWidth',1.5);
yline(0,'r--','s=0');
yline( smc_params.phi,'b:','+\phi');
yline(-smc_params.phi,'b:','-\phi');
xlabel('Temps (s)'); ylabel('s(t)');
title('Surface de glissement'); grid on;

subplot(2,3,5)
plot(t, X(:,4)*180/pi, 'r', 'LineWidth',1.5); hold on;
plot(t, X(:,3)*100,    'b', 'LineWidth',1.5);
legend('\theta\dot (°/s)','\dot{x} (cm/s)','Location','best');
xlabel('Temps (s)'); ylabel('Vitesses');
title('Vitesses'); grid on;

subplot(2,3,6)
angles_test = [5, 10, 15, 20];
colors_ci   = {'b','r','g','m'};
hold on;
for k_ci = 1:length(angles_test)
    X0k = [0; angles_test(k_ci)*pi/180; 0; 0];
    try
        [tk, Xk] = ode45(@(t,Xv) equations_mouvement(t, Xv, ...
            ctrl_mode_glissant(Xv, params, smc_params), params), ...
            [0,8], X0k, options);
        plot(tk, Xk(:,2)*180/pi, colors_ci{k_ci}, 'LineWidth',1.5, ...
            'DisplayName', sprintf('\\theta_0=%d°', angles_test(k_ci)));
    catch
        fprintf('⚠️  Échec θ0=%d°\n', angles_test(k_ci));
    end
end
yline(0,'k--');
xlabel('Temps (s)'); ylabel('\theta (°)');
legend('Location','best');
title('Robustesse - différentes CI'); grid on;

sgtitle('Diagnostic SMC - Pendule Inversé','FontSize',13,'FontWeight','bold');
fprintf('Test terminé ✅\n');