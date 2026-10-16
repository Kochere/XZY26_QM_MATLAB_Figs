%% Fig_S4.m
% Compute QFI for omega, lambda, gamma, and time; derivatives of <H^(1)>;
% and sensitivity quantities B_x(t). Saves results to QFI_results.mat.

clear; clc; close all;

%% Physical parameters and grid
dim = 4;
Is = eye(dim); 
I = eye(2);

sx = [0 1;1 0]; 
sy = [0 -1i;1i 0]; 
sz = [1 0;0 -1];
sm = [0 0;1 0];

om_base = 0.5; 
lam_base = 1.0; 
ga_base = 0.2;
TOL = 1e-10;

d_om = 0.001;
d_lam = 0.001; 
d_ga = 0.001;

down = [0; 1];
psi_0 = kron(down, down);
rho_0 = psi_0 * psi_0';

t = linspace(1e-2, 30+1e-2, 300000);
Nt = length(t);
dt = t(2) - t(1);

% Observable: H^(1) = om*(sz⊗I)
Op = om_base * kron(sz, I);

%% Helper functions (same as original)
function L = build_Lindbladian(om, lam, ga, sx, sz, sm, I, Is, dim)
    H1 = om * kron(sz, I); 
    H2 = om * kron(I, sz); 
    HI = lam * kron(sx, sx);
    Hs = H1 + H2 + HI;
    L_H = -1i * (kron(Is, Hs) - kron(Hs.', Is));
    sm1 = kron(sm, I); 
    sm2 = kron(I, sm);
    L1 = ga * (kron(conj(sm1), sm1) - 0.5 * kron(Is, sm1' * sm1) - 0.5 * kron(sm1.' * conj(sm1), Is));
    L2 = ga * (kron(conj(sm2), sm2) - 0.5 * kron(Is, sm2' * sm2) - 0.5 * kron(sm2.' * conj(sm2), Is));
    L = L_H + L1 + L2;
end

function [R_sort, L_sort, lambda_sort] = sortingEigenvalues(dim, TOL, L)
    dim2 = dim * dim;
    [R, DR] = eig(L);
    [Lmat, DL] = eig(L');
    eig_R = diag(DR);
    eig_L = diag(DL);
    ind_RL = zeros(dim2, 2);
    count = 1;
    for n = 1:dim2
        an = eig_R(n);
        for m = 1:dim2
            bm = eig_L(m);
            if abs(real(an) - real(bm)) < TOL && abs(imag(an) - imag(bm)) < TOL && count <= dim2
                ind_RL(count, :) = [n, m];
                count = count + 1;
            end
        end
    end
    lambda = eig_R(ind_RL(:, 1));
    [~, ind] = sort(real(lambda), 'descend');
    lambda_sort = lambda(ind);
    R_sort = cell(1, dim2);
    L_sort = cell(1, dim2);
    for k = 1:dim2
        Rk = reshape(R(:, ind_RL(ind(k), 1)), dim, dim);
        Lk = reshape(Lmat(:, ind_RL(ind(k), 2)), dim, dim);
        Ck = trace(Lk * Rk);
        if abs(Ck) > 1e-12
            Rk = Rk / sqrt(Ck);
            Lk = Lk / sqrt(Ck);
        end
        R_sort{k} = Rk;
        L_sort{k} = Lk;
    end
end

function F = calc_fidelity(rho, sigma)
    M = rho * sigma;
    lambda = eig(M);
    lambda = real(lambda);
    lambda(lambda < 0) = 0;
    F = abs(sum(sqrt(lambda)))^2;
end

function rho_t = evolve_rho(L, rho_0, t_vec, dim, TOL)
    [R_sort, L_sort, lambda_sort] = sortingEigenvalues(dim, TOL, L);
    Nc = length(lambda_sort);
    c_k = zeros(1, Nc);
    for k = 1:Nc
        c_k(k) = trace(rho_0 * L_sort{k});
    end
    Nt = length(t_vec);
    rho_t = zeros(dim, dim, Nt);
    for n = 1:Nt
        rho = zeros(dim, dim);
        for k = 1:Nc
            rho = rho + c_k(k) * exp(lambda_sort(k) * t_vec(n)) * R_sort{k};
        end
        rho = (rho + rho') / 2;
        rho = rho / trace(rho);
        rho_t(:, :, n) = rho;
    end
end

%% Main computation: QFI for four parameters
fprintf('Computing QFIs...\n');

QFI_om = zeros(1, Nt);
QFI_lam = zeros(1, Nt);
QFI_ga = zeros(1, Nt);
QFI_t = zeros(1, Nt);

% QFI for omega
fprintf('QFI for omega...\n');
L_plus = build_Lindbladian(om_base + d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base - d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
rho_plus = evolve_rho(L_plus, rho_0, t, dim, TOL);
rho_minus = evolve_rho(L_minus, rho_0, t, dim, TOL);
for n = 1:Nt
    F = calc_fidelity(rho_minus(:, :, n), rho_plus(:, :, n));
    if F < 1
        QFI_om(n) = 2 * (1 - F) / (1 + sqrt(F)) / d_om^2;
    end
end

% QFI for lambda
fprintf('QFI for lambda...\n');
L_plus = build_Lindbladian(om_base, lam_base + d_lam, ga_base, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base, lam_base - d_lam, ga_base, sx, sz, sm, I, Is, dim);
rho_plus = evolve_rho(L_plus, rho_0, t, dim, TOL);
rho_minus = evolve_rho(L_minus, rho_0, t, dim, TOL);
for n = 1:Nt
    F = calc_fidelity(rho_minus(:, :, n), rho_plus(:, :, n));
    if F < 1
        QFI_lam(n) = 2 * (1 - F) / (1 + sqrt(F)) / d_lam^2;
    end
end

% QFI for gamma
fprintf('QFI for gamma...\n');
L_plus = build_Lindbladian(om_base, lam_base, ga_base + d_ga, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base, lam_base, ga_base - d_ga, sx, sz, sm, I, Is, dim);
rho_plus = evolve_rho(L_plus, rho_0, t, dim, TOL);
rho_minus = evolve_rho(L_minus, rho_0, t, dim, TOL);
for n = 1:Nt
    F = calc_fidelity(rho_minus(:, :, n), rho_plus(:, :, n));
    if F < 1
        QFI_ga(n) = 2 * (1 - F) / (1 + sqrt(F)) / d_ga^2;
    end
end

% QFI for time (quantum clock)
fprintf('QFI for time...\n');
L_base = build_Lindbladian(om_base, lam_base, ga_base, sx, sz, sm, I, Is, dim);
rho_base = evolve_rho(L_base, rho_0, t, dim, TOL);
for n = 1:Nt
    if n == 1
        F = calc_fidelity(rho_base(:, :, n), rho_base(:, :, n + 1));
        if F < 1
            QFI_t(n) = 8 * (1 - sqrt(F)) / dt^2;
        end
    elseif n == Nt
        F = calc_fidelity(rho_base(:, :, n - 1), rho_base(:, :, n));
        if F < 1
            QFI_t(n) = 8 * (1 - sqrt(F)) / dt^2;
        end
    else
        F = calc_fidelity(rho_base(:, :, n - 1), rho_base(:, :, n + 1));
        if F < 1
            QFI_t(n) = 2 * (1 - F) / (1 + sqrt(F)) / dt^2;
        end
    end
end

%% Observable expectation, variance, and derivatives
fprintf('Computing derivatives and sensitivity...\n');

exp_H = zeros(1, Nt);
exp_H2 = zeros(1, Nt);
for n = 1:Nt
    exp_H(n) = real(trace(Op * rho_base(:, :, n)));
    exp_H2(n) = real(trace(Op * Op * rho_base(:, :, n)));
end
Var_O = exp_H2 - exp_H.^2;

% d<H>/domega
L_plus = build_Lindbladian(om_base + d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base - d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
rho_p = evolve_rho(L_plus, rho_0, t, dim, TOL);
rho_m = evolve_rho(L_minus, rho_0, t, dim, TOL);
exp_p = zeros(1, Nt); exp_m = zeros(1, Nt);
for n = 1:Nt
    exp_p(n) = real(trace(Op * rho_p(:, :, n)));
    exp_m(n) = real(trace(Op * rho_m(:, :, n)));
end
dH_dOm = (exp_p - exp_m) / (2 * d_om);

% d<H>/dlambda
L_plus = build_Lindbladian(om_base, lam_base + d_lam, ga_base, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base, lam_base - d_lam, ga_base, sx, sz, sm, I, Is, dim);
rho_p = evolve_rho(L_plus, rho_0, t, dim, TOL);
rho_m = evolve_rho(L_minus, rho_0, t, dim, TOL);
for n = 1:Nt
    exp_p(n) = real(trace(Op * rho_p(:, :, n)));
    exp_m(n) = real(trace(Op * rho_m(:, :, n)));
end
dH_dLam = (exp_p - exp_m) / (2 * d_lam);

% d<H>/dgamma
L_plus = build_Lindbladian(om_base, lam_base, ga_base + d_ga, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base, lam_base, ga_base - d_ga, sx, sz, sm, I, Is, dim);
rho_p = evolve_rho(L_plus, rho_0, t, dim, TOL);
rho_m = evolve_rho(L_minus, rho_0, t, dim, TOL);
for n = 1:Nt
    exp_p(n) = real(trace(Op * rho_p(:, :, n)));
    exp_m(n) = real(trace(Op * rho_m(:, :, n)));
end
dH_dGa = (exp_p - exp_m) / (2 * d_ga);

% d<H>/dt (numerical)
dH_dt = zeros(1, Nt);
dH_dt(1) = (exp_H(2) - exp_H(1)) / dt;
dH_dt(Nt) = (exp_H(Nt) - exp_H(Nt - 1)) / dt;
for n = 2:Nt - 1
    dH_dt(n) = (exp_H(n + 1) - exp_H(n - 1)) / (2 * dt);
end

%% Sensitivity quantities B_x(t) = (dO/dx)^2 / QFI_x
denom = abs(dH_dLam);
denom(denom < 1e-15) = 1e-15;

ratio_t = abs(dH_dt ./ denom).^2;
ratio_om = abs(dH_dOm ./ denom).^2;
ratio_lam = ones(size(t));
ratio_ga = abs(dH_dGa ./ denom).^2;

min_QFI = 1e-15;
QFI_t_safe = QFI_t;   QFI_t_safe(QFI_t_safe < min_QFI) = inf;
QFI_om_safe = QFI_om; QFI_om_safe(QFI_om_safe < min_QFI) = inf;
QFI_lam_safe = QFI_lam; QFI_lam_safe(QFI_lam_safe < min_QFI) = inf;
QFI_ga_safe = QFI_ga; QFI_ga_safe(QFI_ga_safe < min_QFI) = inf;

quantity_t = ratio_t ./ QFI_t_safe;
quantity_om = ratio_om ./ QFI_om_safe;
quantity_lam = ratio_lam ./ QFI_lam_safe;
quantity_ga = ratio_ga ./ QFI_ga_safe;

quantity_t(~isfinite(quantity_t)) = 1e-20;
quantity_om(~isfinite(quantity_om)) = 1e-20;
quantity_lam(~isfinite(quantity_lam)) = 1e-20;
quantity_ga(~isfinite(quantity_ga)) = 1e-20;

%% Additional: Var(O)/|dH/dLam|^2 for right subplot (black line)
denom_var = abs(dH_dLam).^2;
denom_var(denom_var < 1e-15) = 1e-15;
ratio_var_dLam = Var_O ./ denom_var;
ratio_var_dLam(~isfinite(ratio_var_dLam)) = 1e-20;
ratio_var_dLam(ratio_var_dLam <= 0) = 1e-20;

%% Determine which parameter dominates at each time (for segmented coloring)
all_qty = [quantity_lam; quantity_om; quantity_ga; quantity_t];
[~, max_idx] = max(all_qty);

quantity_lam_max = nan(size(t));
quantity_om_max  = nan(size(t));
quantity_ga_max  = nan(size(t));
quantity_t_max   = nan(size(t));

quantity_lam_max(max_idx == 1) = quantity_lam(max_idx == 1);
quantity_om_max(max_idx == 2)  = quantity_om(max_idx == 2);
quantity_ga_max(max_idx == 3)  = quantity_ga(max_idx == 3);
quantity_t_max(max_idx == 4)   = quantity_t(max_idx == 4);

%% Save results
% save('QFI_results.mat', 't', 'Nt', ...
%      'quantity_lam', 'quantity_om', 'quantity_ga', 'quantity_t', ...
%      'ratio_var_dLam', ...
%      'quantity_lam_max', 'quantity_om_max', 'quantity_ga_max', 'quantity_t_max', ...
%      '-v7.3');

fprintf('Computation done. \n');