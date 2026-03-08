function F = ctrl_mode_glissant(X, params, smc_params)
% CTRL_MODE_GLISSANT  Contrôleur SMC pour le pendule inversé
% Surface basée sur les gains LQR : s = K_lqr * X
% F = F_eq (compense dynamique non linéaire) + F_sw (robustesse)
%
% Entrées :
%   X         : état [x; theta; x_dot; theta_dot]
%   params    : paramètres physiques (M, m, L, I, b, g)
%   smc_params: .K_lqr  gains LQR [1x4]
%               .k      gain robustesse
%               .phi    épaisseur couche limite (anti-chattering)

    % ── Extraction des états ─────────────────────────────────────────────
    x         = X(1);
    theta     = X(2);
    x_dot     = X(3);
    theta_dot = X(4);

    % ── Paramètres physiques ─────────────────────────────────────────────
    M_c = params.M;
    m_p = params.m;
    L   = params.L;
    I   = params.I;
    b   = params.b;
    g   = params.g;

    % ── Paramètres SMC ───────────────────────────────────────────────────
    K_lqr = smc_params.K_lqr;   % vecteur ligne [k1, k2, k3, k4]
    k     = smc_params.k;
    phi   = smc_params.phi;

    % ── Surface de glissement : s = K_lqr * X ───────────────────────────
    s = K_lqr * X;

    % ── Déterminant de la matrice d'inertie ──────────────────────────────
    delta = (M_c + m_p)*(I + m_p*L^2) - (m_p*L*cos(theta))^2;

    % ── Dynamique non linéaire sans F (termes f1, f2) ────────────────────
    % x_ddot     = f1 + g1*F
    % theta_ddot = f2 + g2*F
    f1 = ( (I + m_p*L^2) * (m_p*L*theta_dot^2*sin(theta) - b*x_dot) ...
         + m_p^2*L^2*cos(theta)*g*sin(theta) ) / delta;

    f2 = ( (M_c + m_p)*m_p*g*L*sin(theta) ...
         - m_p*L*cos(theta)*(m_p*L*theta_dot^2*sin(theta) - b*x_dot) ) / delta;

    % ── Coefficients de F dans les accélérations ─────────────────────────
    g1 = (I + m_p*L^2)          / delta;   % dans x_ddot
    g2 = (-m_p*L*cos(theta))    / delta;   % dans theta_ddot

    % ── Dérivée de s : s_dot = A_eq + B_eq*F ────────────────────────────
    % s     = k1*x     + k2*theta     + k3*x_dot     + k4*theta_dot
    % s_dot = k1*x_dot + k2*theta_dot + k3*x_ddot    + k4*theta_ddot
    %       = k1*x_dot + k2*theta_dot + k3*(f1+g1*F) + k4*(f2+g2*F)
    A_eq = K_lqr(1)*x_dot + K_lqr(2)*theta_dot ...
         + K_lqr(3)*f1    + K_lqr(4)*f2;

    B_eq = K_lqr(3)*g1 + K_lqr(4)*g2;

    % Protection division par zéro
    if abs(B_eq) < 1e-8
        B_eq = 1e-8 * sign(B_eq + 1e-20);
    end

    % ── Commande équivalente : impose s_dot = 0 ──────────────────────────
    F_eq = -A_eq / B_eq;

    % ── Terme discontinu avec couche limite (anti-chattering) ────────────
    if abs(s) <= phi
        sat_s = s / phi;        % zone linéaire : pas de discontinuité
    else
        sat_s = sign(s);        % zone glissante : commutation
    end

    F_sw = -(k / abs(B_eq)) * sat_s;

    % ── Commande totale ──────────────────────────────────────────────────
    F = F_eq + F_sw;

    % ── Saturation physique de l'actionneur ──────────────────────────────
    F_max = 20;   % [N]
    F = max(-F_max, min(F_max, F));
end