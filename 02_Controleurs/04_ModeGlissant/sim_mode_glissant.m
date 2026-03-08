% sim_mode_glissant.m
% Simulation du pendule inversé avec contrôle par Mode Glissant (SMC)
% Lancer depuis : 02_Controleurs/04_ModeGlissant/

clear; close all; clc;

% ── 1. Chemins ────────────────────────────────────────────────────────────
base = fileparts(fileparts(mfilename('fullpath')));  % → 02_Controleurs
root = fileparts(base);                               % → PenduleInverse

addpath(fullfile(root, '01_Modelisation'));
addpath(fullfile(base, '04_ModeGlissant'));

% ── 2. Paramètres physiques ───────────────────────────────────────────────
run(fullfile(root, '01_Modelisation', 'parametres_systeme.m'));
params.M = M; params.m = m; params.L = L;
params.I = I; params.b = b; params.g = g;

% ── 3. Gains LQR (surface SMC = surface LQR) ─────────────────────────────
load(fullfile(base, '02_LQR', 'gains_LQR.mat'), 'K');

% ── 4. Paramètres SMC ────────────────────────────────────────────────────
smc_params.K_lqr = K;    % surface s = K*X (même que LQR)
smc_params.k     = 8;    % gain robustesse (terme discontinu)
smc_params.phi   = 0.3;  % couche limite anti-chattering

% ── 5. Condition initiale ────────────────────────────────────────────────
X0 = [0; 10*pi/180; 0; 0];   % 10° d'inclinaison initiale

% ── 6. Simulation ODE45 ──────────────────────────────────────────────────
t_span  = [0, 10];
options = odeset('RelTol',1e-6, 'AbsTol',1e-8, 'MaxStep', 0.01);

[t, X] = ode45(@(t,Xv) equations_mouvement(t, Xv, ...
    ctrl_mode_glissant(Xv, params, smc_params), params), ...
    t_span, X0, options);

% ── 7. Recalcul commande et surface ──────────────────────────────────────
F_hist = zeros(length(t), 1);
s_hist = zeros(length(t), 1);

for i = 1:length(t)
    F_hist(i) = ctrl_mode_glissant(X(i,:)', params, smc_params);
    s_hist(i) = smc_params.K_lqr * X(i,:)';
end

% ── 8. Métriques de performance ───────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('   RÉSULTATS SMC - Pendule Inversé\n');
fprintf('========================================\n');

idx1 = find(abs(X(:,2)) < 1*pi/180, 1, 'first');
if ~isempty(idx1)
    fprintf('✅ Stabilisation θ < 1°   : %.2f s\n', t(idx1));
else
    fprintf('⚠️  θ non stabilisé à 1° sur 10s\n');
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
fprintf('========================================\n\n');

% ── 9. Figures ────────────────────────────────────────────────────────────
figure('Name','SMC - Pendule Inversé', ...
       'NumberTitle','off', ...
       'Position',[50 50 1300 700]);

% Angle du pendule
subplot(2,3,1)
plot(t, X(:,2)*180/pi, 'r', 'LineWidth', 2);
yline(0,  'k--', 'Réf');
yline( 1, 'b:',  '+1°');
yline(-1, 'b:',  '-1°');
xlabel('Temps (s)'); ylabel('\theta (°)');
title('Angle du pendule'); grid on;

% Position du chariot
subplot(2,3,2)
plot(t, X(:,1)*100, 'b', 'LineWidth', 2);
yline(0, 'k--', 'Réf');
xlabel('Temps (s)'); ylabel('x (cm)');
title('Position du chariot'); grid on;

% Commande
subplot(2,3,3)
plot(t, F_hist, 'm', 'LineWidth', 1.5);
yline( 20, 'r:', '+F_{max}');
yline(-20, 'r:', '-F_{max}');
xlabel('Temps (s)'); ylabel('F (N)');
title('Force de commande'); grid on;

% Surface de glissement
subplot(2,3,4)
plot(t, s_hist, 'k', 'LineWidth', 1.5);
yline(0,                  'r--', 's=0');
yline( smc_params.phi,    'b:',  '+\phi');
yline(-smc_params.phi,    'b:',  '-\phi');
xlabel('Temps (s)'); ylabel('s(t)');
title('Surface de glissement'); grid on;

% Vitesses
subplot(2,3,5)
plot(t, X(:,4)*180/pi, 'r', 'LineWidth', 1.5); hold on;
plot(t, X(:,3)*100,    'b', 'LineWidth', 1.5);
xlabel('Temps (s)'); ylabel('Vitesses');
legend('\theta\dot (°/s)', '\dot{x} (cm/s)', 'Location','best');
title('Vitesses'); grid on;

% Robustesse : différentes conditions initiales
subplot(2,3,6)
angles_test = [5, 10, 15, 20];
colors_ci   = {'b','r','g','m'};
hold on;
for k_ci = 1:length(angles_test)
    X0k = [0; angles_test(k_ci)*pi/180; 0; 0];
    try
        [tk, Xk] = ode45(@(t,Xv) equations_mouvement(t, Xv, ...
            ctrl_mode_glissant(Xv, params, smc_params), params), ...
            [0, 8], X0k, options);
        plot(tk, Xk(:,2)*180/pi, colors_ci{k_ci}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('\\theta_0=%d°', angles_test(k_ci)));
    catch
        fprintf('⚠️  Échec simulation θ0=%d°\n', angles_test(k_ci));
    end
end
yline(0, 'k--');
xlabel('Temps (s)'); ylabel('\theta (°)');
legend('Location','best');
title('Robustesse - différentes CI'); grid on;

sgtitle('Contrôle par Mode Glissant (SMC)', ...
        'FontSize', 14, 'FontWeight', 'bold');