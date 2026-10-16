%% Fig_2.m
% Compute for two initial states: estimator variance, our bound, QCRB
% Results saved in workspace: t, (Var_A, Bound_A, QCRB_A), (Var_B, Bound_B, QCRB_B)

clear; close all;

%% Parameters
sx = [0 1; 1 0];
sz = [1 0; 0 -1];
omega0 = 0.5;
eta0   = 1;

t = linspace(1e-2, 20+1e-2, 6000);
Nt = length(t);
dt = t(2) - t(1);
d_theta = 0.001;

O = sz;

%% Initial states
% State A: pure ground state |g> = [0;1]
% State B: mixed state with coherence
psi0_A = [0; 1];
rho0_A = psi0_A * psi0_A';

rho0_B = [0.6, 0.4+0.2i; 0.4-0.2i, 0.4];

%% Compute for both states
fprintf('Computing for initial state A ...\n');
[Var_A, Bound_A, QCRB_A] = compute_for_state(rho0_A, omega0, eta0, t, dt, d_theta, sx, sz, O);

fprintf('Computing for initial state B ...\n');
[Var_B, Bound_B, QCRB_B] = compute_for_state(rho0_B, omega0, eta0, t, dt, d_theta, sx, sz, O);

fprintf('Computation done. \n');

%% Function definitions
function [Var_eta, Bound_eta, QCRB_eta] = compute_for_state(rho0, omega0, eta0, t, dt, d_theta, sx, sz, O)
    % Compute variance, bound, and QCRB for a given initial state
    Nt = length(t);
    
    % <O> and Var(O)
    [Ovals, varO] = compute_mean_and_var(omega0, eta0, rho0, t, O, sx, sz);
    
    % QFI for omega (via fidelity)
    rho_m_om = evolve_state(omega0-d_theta, eta0, rho0, t, sx, sz);
    rho_p_om = evolve_state(omega0+d_theta, eta0, rho0, t, sx, sz);
    QFI_omega = zeros(1, Nt);
    for ii = 1:Nt
        F = fidelity(rho_m_om(:,:,ii), rho_p_om(:,:,ii));
        QFI_omega(ii) = 2 * (1 - sqrt(F)) / d_theta^2;
    end
    QFI_omega = max(QFI_omega, 0);
    
    % QFI for eta
    rho_m_eta = evolve_state(omega0, eta0-d_theta, rho0, t, sx, sz);
    rho_p_eta = evolve_state(omega0, eta0+d_theta, rho0, t, sx, sz);
    QFI_eta = zeros(1, Nt);
    for ii = 1:Nt
        F = fidelity(rho_m_eta(:,:,ii), rho_p_eta(:,:,ii));
        QFI_eta(ii) = 2 * (1 - sqrt(F)) / d_theta^2;
    end
    QFI_eta = max(QFI_eta, 0);
    QCRB_eta = 1 ./ QFI_eta;
    
    % QFI for time (needed for the bound)
    rho_base = evolve_state(omega0, eta0, rho0, t, sx, sz);
    QFI_t = zeros(1, Nt);
    for ii = 1:Nt
        if ii == 1
            F = fidelity(rho_base(:,:,ii), rho_base(:,:,ii+1));
            QFI_t(ii) = 8 * (1 - sqrt(F)) / dt^2;
        elseif ii == Nt
            F = fidelity(rho_base(:,:,ii-1), rho_base(:,:,ii));
            QFI_t(ii) = 8 * (1 - sqrt(F)) / dt^2;
        else
            F = fidelity(rho_base(:,:,ii-1), rho_base(:,:,ii+1));
            QFI_t(ii) = 2 * (1 - F) / (1 + sqrt(F)) / dt^2;
        end
    end
    QFI_t = max(QFI_t, 0);
    
    % Derivative d<O>/deta
    O_p = compute_expectation(omega0, eta0+d_theta, rho0, t, O, sx, sz);
    O_m = compute_expectation(omega0, eta0-d_theta, rho0, t, O, sx, sz);
    dO_deta = (O_p - O_m) / (2*d_theta);
    Var_eta = varO ./ (dO_deta.^2 + 1e-12);   % error propagation
    
    % Derivative d<O>/domega
    O_p_om = compute_expectation(omega0+d_theta, eta0, rho0, t, O, sx, sz);
    O_m_om = compute_expectation(omega0-d_theta, eta0, rho0, t, O, sx, sz);
    dO_domega = (O_p_om - O_m_om) / (2*d_theta);
    
    % Derivative d<O>/dt
    dO_dt = zeros(1, Nt);
    dO_dt(1) = (Ovals(2) - Ovals(1)) / dt;
    dO_dt(end) = (Ovals(end) - Ovals(end-1)) / dt;
    dO_dt(2:end-1) = (Ovals(3:end) - Ovals(1:end-2)) / (2*dt);
    
    % Bound: maximum of three terms
    denom = abs(dO_deta);
    denom(denom < 1e-15) = 1e-15;
    ratio_eta = ones(size(t));                      
    ratio_om  = abs(dO_domega ./ denom).^2;
    ratio_t   = abs(dO_dt ./ denom).^2;
    
    % Safeguard against zero QFI
    minQ = 1e-15;
    QFI_om_safe = QFI_omega;  QFI_om_safe(QFI_om_safe < minQ) = inf;
    QFI_eta_safe = QFI_eta;   QFI_eta_safe(QFI_eta_safe < minQ) = inf;
    QFI_t_safe   = QFI_t;     QFI_t_safe(QFI_t_safe < minQ)     = inf;
    
    quantity_eta = ratio_eta ./ QFI_eta_safe;
    quantity_om  = ratio_om  ./ QFI_om_safe;
    quantity_t   = ratio_t   ./ QFI_t_safe;
    
    Bound_eta = max([quantity_eta; quantity_om; quantity_t], [], 1);
    Bound_eta(~isfinite(Bound_eta)) = 1e-20;   % replace non-finite values
end

% -------------------------------------------------------------------------
% Helper functions
% -------------------------------------------------------------------------

function rho_t = evolve_state(omega, eta, rho0, t_vec, sx, sz)
    % Evolve state under H = omega*sz + eta*sx
    H = omega * sz + eta * sx;
    Nt = length(t_vec);
    rho_t = zeros(2, 2, Nt);
    for ii = 1:Nt
        U = expm(-1i * H * t_vec(ii));
        rho_t(:,:,ii) = U * rho0 * U';
    end
end

function Ovals = compute_expectation(omega, eta, rho0, t_vec, O, sx, sz)
    % Expectation value <O>(t) = Tr(O rho(t))
    H = omega * sz + eta * sx;
    Nt = length(t_vec);
    Ovals = zeros(1, Nt);
    for ii = 1:Nt
        U = expm(-1i * H * t_vec(ii));
        rho = U * rho0 * U';
        Ovals(ii) = real(trace(O * rho));
    end
end

function [Ovals, varO] = compute_mean_and_var(omega, eta, rho0, t_vec, O, sx, sz)
    % Compute <O> and variance Var(O) = <O^2> - <O>^2
    H = omega * sz + eta * sx;
    Nt = length(t_vec);
    Ovals = zeros(1, Nt);
    varO = zeros(1, Nt);
    for ii = 1:Nt
        U = expm(-1i * H * t_vec(ii));
        rho = U * rho0 * U';
        Ovals(ii) = real(trace(O * rho));
        O2 = real(trace(O * O * rho));
        varO(ii) = max(O2 - Ovals(ii)^2, 0);
    end
end

function F = fidelity(rho, sigma)
    % Fidelity F(rho, sigma) = [Tr(sqrt(sqrt(rho)*sigma*sqrt(rho)))]^2
    M = rho * sigma;
    lam = eig(M);
    lam = real(lam);
    lam(lam < 0) = 0;
    F = (sum(sqrt(lam)))^2;
    F = min(real(F), 1);
end