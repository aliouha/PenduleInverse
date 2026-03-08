function F = ctrl_backstepping(X, params, bs_params)
% CTRL_BACKSTEPPING - Pendule inversé
% Principe : annulation de dynamique non linéaire + retour d'état
% Surface : s = K_lqr * X  (même que SMC, prouvé stable)
% F = F_eq (annule dynamique) + F_amort (amortissement)

    x         = X(1);
    theta     = X(2);
    x_dot     = X(3);
    theta_dot = X(4);

    M_c = params.M;  m_p = params.m;
    L   = params.L;  I   = params.I;
    b   = params.b;  g   = params.g;

    K   = bs_params.K_lqr;   % gains LQR [1x4]
    kb  = bs_params.kb;      % gain backstepping (amortissement)

    % ── Déterminant ────────────────────────────────────────────────────────
    delta = (M_c + m_p)*(I + m_p*L^2) - (m_p*L*cos(theta))^2;

    % ── Dynamique non linéaire : accél = f + g*F ───────────────────────────
    f1 = ( (I + m_p*L^2)*(m_p*L*theta_dot^2*sin(theta) - b*x_dot) ...
         + m_p^2*L^2*cos(theta)*g*sin(theta) ) / delta;

    f2 = ( (M_c + m_p)*m_p*g*L*sin(theta) ...
         - m_p*L*cos(theta)*(m_p*L*theta_dot^2*sin(theta) - b*x_dot) ) / delta;

    g1 = (I + m_p*L^2)       / delta;
    g2 = (-m_p*L*cos(theta)) / delta;

    % ── Surface s = K*X ────────────────────────────────────────────────────
    s = K * X;

    % ── Dérivée de s : s_dot = A + B*F ─────────────────────────────────────
    A = K(1)*x_dot + K(2)*theta_dot + K(3)*f1 + K(4)*f2;
    B = K(3)*g1    + K(4)*g2;

    if abs(B) < 1e-8, B = 1e-8*sign(B + 1e-20); end

    % ── Commande : annule dynamique + amortit la surface ───────────────────
    % On veut s_dot = -kb * s  (convergence exponentielle de s vers 0)
    % A + B*F = -kb*s  →  F = (-kb*s - A) / B
    F = (-kb*s - A) / B;

    % ── Saturation ─────────────────────────────────────────────────────────
    F_max = 20;
    F = max(-F_max, min(F_max, F));
end