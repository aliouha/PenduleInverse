% comparaison_tous_controleurs.m
% Comparaison finale : PID vs LQR vs LQG-Kalman vs SMC vs Backstepping
% PID testé à 2° (limite), autres à 10°
% Lancer depuis : 03_Analyse_Comparaison/

clear; close all; clc;

fprintf('========================================\n');
fprintf('   COMPARAISON TOUS CONTROLEURS\n');
fprintf('========================================\n\n');

% ── Chemins ───────────────────────────────────────────────────────────────
root = 'E:\Mes-CD\Matlab-Projects\PenduleInverse\PenduleInverse';
ctrl = fullfile(root, '02_Controleurs');

addpath(fullfile(root, '01_Modelisation'));
addpath(fullfile(ctrl, '04_ModeGlissant'));
addpath(fullfile(ctrl, '05_Backstepping'));

% ── Paramètres physiques ──────────────────────────────────────────────────
run(fullfile(root, '01_Modelisation', 'parametres_systeme.m'));
params = struct('M',M,'m',m,'L',L,'I',I,'b',b,'g',g);
load(fullfile(root, '01_Modelisation', 'modele_lineaire.mat'), 'A_num', 'B_num');

% ── Chargement des gains ──────────────────────────────────────────────────
load(fullfile(ctrl, '01_PID',        'gains_PID.mat'));
load(fullfile(ctrl, '02_LQR',        'gains_LQR.mat'), 'K');
load(fullfile(ctrl, '06_LQG_kalman', 'gains_LQG_Kalman.mat'), 'L_kalman');

gains_pid = struct('Kp_theta',Kp_theta,'Kd_theta',Kd_theta,...
                   'Ki_theta',Ki_theta,'Kp_x',Kp_x,'Kd_x',Kd_x,'F_max',F_max);

smc_params.K_lqr = K; smc_params.k = 8; smc_params.phi = 0.3;
bs_params.K_lqr  = K; bs_params.kb  = 5;

C_obs = [1 0 0 0; 0 1 0 0];
Rk    = diag([0.001, 0.001]);

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.01);

% ── Conditions initiales ──────────────────────────────────────────────────
X0_pid    = [0; 2*pi/180;  0; 0];    % PID : 2° (limite du PID)
X0_others = [0; 10*pi/180; 0; 0];    % Autres : 10°

% ── 1. PID à 2° ───────────────────────────────────────────────────────────
fprintf('Simulation PID (2 deg)...  ');
X0_pid_ext = [X0_pid; 0];
[t_pid, X_pid_ext] = ode45(@(t,Xe) pendule_PID_ext(t,Xe,params,gains_pid), ...
    [0 10], X0_pid_ext, opts);
X_pid = X_pid_ext(:,1:4);
F_pid = zeros(length(t_pid),1);
for i=1:length(t_pid)
    F_pid(i) = commande_PID_ext(X_pid_ext(i,:)', gains_pid);
end
pid_ok = max(abs(X_pid(:,2))) < 90*pi/180;
fprintf('OK - %s\n', ternaire(pid_ok,'STABLE','INSTABLE'));

% ── 2. LQR à 10° ─────────────────────────────────────────────────────────
fprintf('Simulation LQR (10 deg)... ');
[t_lqr, X_lqr] = ode45(@(t,X) equations_mouvement(t,X,-K*X,params), ...
    [0 10], X0_others, opts);
F_lqr = (-K * X_lqr')';
fprintf('OK\n');

% ── 3. LQG-Kalman à 10° ──────────────────────────────────────────────────
fprintf('Simulation LQG-Kalman...   ');
dt_lqg=0.005; t_lqg=(0:dt_lqg:10)'; N_lqg=length(t_lqg);
Xr_lqg=zeros(N_lqg,4); Xe_lqg=zeros(N_lqg,4); F_lqg=zeros(N_lqg,1);
Xr_lqg(1,:)=X0_others'; Xe_lqg(1,:)=X0_others';
sx=sqrt(Rk(1,1)); st=sqrt(Rk(2,2));
for i=1:N_lqg-1
    F_lqg(i)=-K*Xe_lqg(i,:)';
    y=C_obs*Xr_lqg(i,:)'+[sx*randn;st*randn];
    dXr=equations_mouvement(t_lqg(i),Xr_lqg(i,:)',F_lqg(i),params);
    Xr_lqg(i+1,:)=Xr_lqg(i,:)+dXr'*dt_lqg;
    dXe=A_num*Xe_lqg(i,:)'+B_num*F_lqg(i)+L_kalman*(y-C_obs*Xe_lqg(i,:)');
    Xe_lqg(i+1,:)=Xe_lqg(i,:)+dXe'*dt_lqg;
end
fprintf('OK\n');

% ── 4. SMC à 10° ─────────────────────────────────────────────────────────
fprintf('Simulation SMC (10 deg)... ');
[t_smc, X_smc] = ode45(@(t,X) equations_mouvement(t,X,...
    ctrl_mode_glissant(X,params,smc_params),params), [0 10], X0_others, opts);
F_smc = arrayfun(@(i) ctrl_mode_glissant(X_smc(i,:)',params,smc_params), ...
    1:length(t_smc))';
fprintf('OK\n');

% ── 5. Backstepping à 10° ────────────────────────────────────────────────
fprintf('Simulation Backstepping... ');
[t_bs, X_bs] = ode45(@(t,X) equations_mouvement(t,X,...
    ctrl_backstepping(X,params,bs_params),params), [0 10], X0_others, opts);
F_bs = arrayfun(@(i) ctrl_backstepping(X_bs(i,:)',params,bs_params), ...
    1:length(t_bs))';
fprintf('OK\n\n');

% ── Métriques ─────────────────────────────────────────────────────────────
noms       = {'PID','LQR','LQG-Kalman','SMC','Backstepping'};
theta0_str = {'2 deg','10 deg','10 deg','10 deg','10 deg'};
temps_l    = {t_pid, t_lqr, t_lqg, t_smc, t_bs};
theta_l    = {X_pid(:,2), X_lqr(:,2), Xr_lqg(:,2), X_smc(:,2), X_bs(:,2)};
x_l        = {X_pid(:,1), X_lqr(:,1), Xr_lqg(:,1), X_smc(:,1), X_bs(:,1)};
F_l        = {F_pid, F_lqr, F_lqg, F_smc, F_bs};

fprintf('%-14s %-8s %-12s %-12s %-10s %-10s\n',...
    'Controleur','theta0','t_stab(s)','x_max(cm)','F_max(N)','Statut');
fprintf('%s\n', repmat('-',1,72));

t_stab_v=[]; x_max_v=[]; F_max_v=[];
for k=1:5
    th=theta_l{k}; tt=temps_l{k};
    idx=find(abs(th)<1*pi/180,1,'first');
    ts=NaN; if ~isempty(idx), ts=tt(idx); end
    xm=max(abs(x_l{k}))*100;
    fm=max(abs(F_l{k}));
    stable = max(abs(th))<90*pi/180;
    if stable, statut='STABLE'; else, statut='INSTABLE'; end
    fprintf('%-14s %-8s %-12s %-12.1f %-10.1f %-10s\n',...
        noms{k}, theta0_str{k}, ...
        ternaire(~isnan(ts),sprintf('%.2f',ts),'N/A'), xm, fm, statut);
    t_stab_v(k)=ts; x_max_v(k)=xm; F_max_v(k)=fm;
end
fprintf('\n* PID teste a 2 deg (instable au-dela)\n\n');

% ── Couleurs et styles ────────────────────────────────────────────────────
C = {[0 0.45 0.74],[0.47 0.67 0.19],[0.93 0.69 0.13],[0.85 0.33 0.1],[0.49 0.18 0.56]};
S = {'-','--',':','-.','--'};
lw = 2;

figure('Name','Comparaison Controleurs','NumberTitle','off',...
       'Position',[20 20 1500 850]);

% ── Subplot 1 : Angle theta (PID séparé) ─────────────────────────────────
subplot(2,3,1); hold on;
% PID sur axe secondaire car diverge
yyaxis left
plot(t_pid, X_pid(:,2)*180/pi, S{1},'Color',C{1},'LineWidth',lw);
ylabel('theta PID (deg)','Color',C{1});
yyaxis right
plot(t_lqr, X_lqr(:,2)*180/pi, S{2},'Color',C{2},'LineWidth',lw);
plot(t_lqg, Xr_lqg(:,2)*180/pi,S{3},'Color',C{3},'LineWidth',lw);
plot(t_smc, X_smc(:,2)*180/pi, S{4},'Color',C{4},'LineWidth',lw);
plot(t_bs,  X_bs(:,2)*180/pi,  S{5},'Color',C{5},'LineWidth',lw);
ylabel('theta LQR/LQG/SMC/BS (deg)');
yline(0,'k--','LineWidth',1);
xlabel('Temps (s)');
title('Angle theta (PID=2 deg, autres=10 deg)'); grid on;
legend(noms,'Location','northeast','FontSize',7);

% ── Subplot 2 : Angle des contrôleurs stables seulement ──────────────────
subplot(2,3,2); hold on;
plot(t_lqr, X_lqr(:,2)*180/pi, S{2},'Color',C{2},'LineWidth',lw);
plot(t_lqg, Xr_lqg(:,2)*180/pi,S{3},'Color',C{3},'LineWidth',lw);
plot(t_smc, X_smc(:,2)*180/pi, S{4},'Color',C{4},'LineWidth',lw);
plot(t_bs,  X_bs(:,2)*180/pi,  S{5},'Color',C{5},'LineWidth',lw);
yline(0,'k--'); yline(1,'k:'); yline(-1,'k:');
xlabel('Temps (s)'); ylabel('theta (deg)');
title('Angle - LQR/LQG/SMC/BS (theta0=10 deg)'); grid on;
legend(noms(2:5),'Location','northeast','FontSize',8);

% ── Subplot 3 : Position x ───────────────────────────────────────────────
subplot(2,3,3); hold on;
plot(t_lqr, X_lqr(:,1)*100, S{2},'Color',C{2},'LineWidth',lw);
plot(t_lqg, Xr_lqg(:,1)*100,S{3},'Color',C{3},'LineWidth',lw);
plot(t_smc, X_smc(:,1)*100, S{4},'Color',C{4},'LineWidth',lw);
plot(t_bs,  X_bs(:,1)*100,  S{5},'Color',C{5},'LineWidth',lw);
yline(0,'k--');
xlabel('Temps (s)'); ylabel('x (cm)');
title('Position chariot (theta0=10 deg)'); grid on;
legend(noms(2:5),'Location','northeast','FontSize',8);

% ── Subplot 4 : Commande ─────────────────────────────────────────────────
subplot(2,3,4); hold on;
plot(t_pid, F_pid, S{1},'Color',C{1},'LineWidth',lw);
plot(t_lqr, F_lqr, S{2},'Color',C{2},'LineWidth',lw);
plot(t_lqg, F_lqg, S{3},'Color',C{3},'LineWidth',lw);
plot(t_smc, F_smc, S{4},'Color',C{4},'LineWidth',lw);
plot(t_bs,  F_bs,  S{5},'Color',C{5},'LineWidth',lw);
xlabel('Temps (s)'); ylabel('F (N)');
title('Commande (force)'); grid on;
legend(noms,'Location','northeast','FontSize',7);

% ── Subplot 5 : Barplot temps stabilisation ───────────────────────────────
subplot(2,3,5);
valid_idx = find(~isnan(t_stab_v));
vals = t_stab_v; vals(isnan(vals)) = 0;
b = bar(1:5, vals, 0.6);
b.FaceColor = 'flat';
for k=1:5
    if ~isnan(t_stab_v(k))
        b.CData(k,:) = C{k};
    else
        b.CData(k,:) = [0.85 0.85 0.85];
    end
end
set(gca,'XTick',1:5,'XTickLabel',noms,'XTickLabelRotation',20,'FontSize',8);
ylabel('Temps (s)'); title('Temps stabilisation theta < 1 deg'); grid on;
% Annotation PID
text(1, 0.05, 'Non', 'HorizontalAlignment','center','FontSize',8,'Color','r');
text(1, 0.02, 'stab.','HorizontalAlignment','center','FontSize',8,'Color','r');

% ── Subplot 6 : Barplot force max ─────────────────────────────────────────
subplot(2,3,6);
b2 = bar(1:5, F_max_v, 0.6);
b2.FaceColor = 'flat';
for k=1:5, b2.CData(k,:) = C{k}; end
set(gca,'XTick',1:5,'XTickLabel',noms,'XTickLabelRotation',20,'FontSize',8);
ylabel('F max (N)'); title('Force maximale appliquee'); grid on;

sgtitle({'Comparaison des Controleurs - Pendule Inverse',...
         'PID: theta0=2 deg  |  LQR, LQG, SMC, Backstepping: theta0=10 deg'},...
        'FontSize',12,'FontWeight','bold');

% ── Fonctions locales ──────────────────────────────────────────────────────
function F = commande_PID_ext(Xe, gains)
    theta=Xe(2); x=Xe(1); x_dot=Xe(3); theta_dot=Xe(4); int_th=Xe(5);
    F = gains.Kp_theta*(0-theta) + gains.Kd_theta*(0-theta_dot) ...
      + gains.Ki_theta*int_th ...
      + gains.Kp_x*(0-x) + gains.Kd_x*(0-x_dot);
    F = max(min(F, gains.F_max), -gains.F_max);
end

function dXe = pendule_PID_ext(t, Xe, params, gains)
    F  = commande_PID_ext(Xe, gains);
    dX = equations_mouvement(t, Xe(1:4), F, params);
    dXe = [dX; Xe(2)];
end

function r = ternaire(cond, a, b)
    if cond, r=a; else, r=b; end
end