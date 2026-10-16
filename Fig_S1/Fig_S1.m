%% Fig_S1.m
% Compute QFI for temperature estimation in a dissipative qubit.
% Saves variables to workspace: t, QFI_true, QFI_num, F_ss_theory.

clear; clc; close all;

%% Parameters (matching the original)
T_base = 0.5;           % sample temperature
omega = 1.0;            % level splitting
gamma = 1.0;            % coupling to bath
gamma0 = 0.2;           % pure dephasing rate

TOL = 1e-12;
d_T = 0.001;            % finite-difference step for temperature

% Initial conditions:
p0 = -0.2;
f0 = 0.4;

% Time grid (avoid t=0)
t_max = 8;
Nt = 50;
t = linspace(0.01, t_max, Nt);

%% Derived physical parameters
beta = 1/T_base;
N_th = 1/(exp(omega/T_base) - 1);
gamma_plus = gamma * N_th;
gamma_minus = gamma * (N_th + 1);
gamma_p = gamma_plus + gamma_minus;
gamma_m = gamma_plus - gamma_minus;

% Steady-state values
p_ss = gamma_m / (2*gamma_p);
mu_plus = 0.5 + p_ss;
mu_minus = 0.5 - p_ss;

% Time-dependent exponentials
exp_gammat = exp(-gamma_p * t);
exp_2gammat = exp(-2*gamma_p * t);
exp_4gammat = exp(-(4*gamma0 + gamma_p)*t);

%% Analytical solution for QFI
p_t = p_ss + (p0 - p_ss)*exp_gammat;
f_t = f0 * exp(-1i*omega*t) .* exp(-(2*gamma0 + gamma_p/2)*t);
f_t_abs_sq = abs(f_t).^2;

H_val = (exp_gammat + gamma_p.*t - 1)./(gamma_p.*t);

term1_code = abs(2*p_ss*H_val - p_t).^2 .* f_t_abs_sq;
numerator_term2 = (p_t.*(p_ss*H_val - p_t) - 0.5*f_t_abs_sq).^2;
term2_code = numerator_term2 / (mu_plus * mu_minus);

denom = gamma_p^2 * T_base^4 * p_ss^2 * (1 - 4*mu_plus*mu_minus);
prefactor_code = (4 * omega^2 * t.^2 * gamma_plus^2 * gamma_minus^2) / denom;
QFI_true = prefactor_code .* (term1_code + term2_code);

%% Numerical QFI via fidelity (finite difference in temperature)
fprintf('Computing numerical QFI via fidelity...\n');

% Initial density matrix: rho = 0.5*I + p0*sz + f0*sx
rho_0 = [0.5+p0, f0; f0, 0.5-p0];

% Build Lindbladian at base temperature
L_base = build_L(omega, gamma, gamma0, T_base);

% Evolve at base temperature for all times
rho_base = cell(1, Nt);
for n = 1:Nt
    rho_base{n} = evolve_rho(L_base, rho_0, t(n), TOL);
end

% Perturbed Lindbladians
L_plus = build_L(omega, gamma, gamma0, T_base + d_T);
L_minus = build_L(omega, gamma, gamma0, T_base - d_T);

QFI_num = zeros(1, Nt);
for n = 1:Nt
    rho_p = evolve_rho(L_plus, rho_0, t(n), TOL);
    rho_m = evolve_rho(L_minus, rho_0, t(n), TOL);
    F = fidelity(rho_m, rho_p);
    F = min(F, 1.0);
    % Symmetric finite-difference formula for QFI
    QFI_num(n) = 2*(1 - F) / (1 + sqrt(F)) / d_T^2;
end

% Steady-state QFI (theoretical)
F_ss_theory = omega^2/(4*T_base^4*cosh(omega/(2*T_base))^2);
% Alternative (unused in plot): F_ss = (omega^2*gamma_minus*gamma_plus)/(T_base^4*gamma_p^2);

%% Save variables to workspace (used by Plot_Fig_S1)
fprintf('Computation done. \n');

%% Helper functions
function L = build_L(omega, gamma, gamma0, T)
    % Build Lindbladian superoperator (4x4) for the qubit.
    if T > 0
        N = 1/(exp(omega/T) - 1);
    else
        N = 0;
    end
    
    gamma_p_loc = gamma * N;
    gamma_m_loc = gamma * (N+1);
    
    sz = [1 0; 0 -1];
    sp = [0 1; 0 0];
    sm = [0 0; 1 0];
    H = 0.5*omega*sz;
    I2 = eye(2);
    
    L_H = -1i*(kron(I2, H) - kron(H.', I2));
    L1 = gamma_p_loc*(kron(conj(sp), sp) - 0.5*kron(I2, sp'*sp) - 0.5*kron(sp.'*conj(sp), I2));
    L2 = gamma_m_loc*(kron(conj(sm), sm) - 0.5*kron(I2, sm'*sm) - 0.5*kron(sm.'*conj(sm), I2));
    L3 = gamma0*(kron(conj(sz), sz) - 0.5*kron(I2, sz'*sz) - 0.5*kron(sz.'*conj(sz), I2));
    
    L = L_H + L1 + L2 + L3;
end

function rho_t = evolve_rho(L, rho_0, t_vec, TOL)
    % Evolve density matrix under Lindbladian L for time vector t_vec.
    dim = 2; dim2 = 4;
    [R, DR] = eig(L);
    [Lmat, DL] = eig(L');
    eig_R = diag(DR); eig_L = diag(DL);
    
    % Match eigenvalues between L and L' (for biorthogonal basis)
    ind_RL = zeros(dim2, 2); count = 1;
    for n = 1:dim2
        for m = 1:dim2
            if abs(real(eig_R(n))-real(eig_L(m)))<TOL && abs(imag(eig_R(n))-imag(eig_L(m)))<TOL && count<=dim2
                ind_RL(count, :) = [n, m]; count = count + 1;
            end
        end
    end
    
    lambda = eig_R(ind_RL(:,1));
    [~, ind] = sort(real(lambda), 'descend');
    lambda_sort = lambda(ind);
    
    % Build eigenmatrices
    R_mat = cell(1, dim2); L_mat = cell(1, dim2);
    for k = 1:dim2
        Rk = reshape(R(:,ind_RL(ind(k),1)), dim, dim);
        Lk = reshape(Lmat(:,ind_RL(ind(k),2)), dim, dim);
        Ck = trace(Lk * Rk);
        if abs(Ck) > 1e-12
            Rk = Rk/sqrt(Ck);
            Lk = Lk/sqrt(Ck);
        end
        R_mat{k} = Rk;
        L_mat{k} = Lk;
    end
    
    % Expansion coefficients
    c_k = zeros(1, dim2);
    for k = 1:dim2
        c_k(k) = trace(rho_0 * L_mat{k});
    end
    
    % Time evolution
    Nt = length(t_vec);
    rho_t = zeros(2, 2, Nt);
    for n = 1:Nt
        rho = zeros(2, 2);
        for k = 1:dim2
            rho = rho + c_k(k) * exp(lambda_sort(k)*t_vec(n)) * R_mat{k};
        end
        rho = (rho + rho')/2;
        rho = rho / trace(rho);
        rho_t(:,:,n) = rho;
    end
end

function F = fidelity(rho, sigma)
    % Fidelity for 2x2 density matrices.
    M = rho * sigma;
    lambda = real(eig(M));
    lambda(lambda < 0) = 0;
    F = abs(sum(sqrt(lambda)))^2;
end