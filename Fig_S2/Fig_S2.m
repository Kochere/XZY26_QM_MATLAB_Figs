%% Fig_S2.m
% Compute QFI(t) for time, omega, lambda, gamma in an open two-spin system.
% Saves variables: t, QFI_t, QFI_om, QFI_lam, QFI_ga.

clc; clear; close all;

%% Basic settings
dim = 4;
Is = eye(dim); I = eye(2);
% Pauli operators
sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];
sp = 0.5*(sx+1i*sy); sm = 0.5*(sx-1i*sy);

% Parameters
om_base = 0.5; 
lam_base = 1; 
ga_base = 0.2;
TOL = 1e-10;

% Initial state |↓↓>
down = [0; 1];
psi_0 = kron(down, down);
rho_0 = psi_0 * psi_0';

% Time grid
t = linspace(0, 30, 30000);
Nt = length(t);

%% Compute QFI for four parameters via fidelity finite difference
d_om = 0.001;
d_lam = 0.001;
d_ga = 0.001;
dt = t(2) - t(1);

% Preallocate
QFI_om = zeros(1, Nt);
QFI_lam = zeros(1, Nt);
QFI_ga = zeros(1, Nt);
QFI_t = zeros(1, Nt);

% 1. QFI for omega
fprintf('Calculating QFI for omega...\n');
L_plus = build_Lindbladian(om_base+d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base-d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
rho_plus = evolve_rho_t(L_plus, rho_0, t, dim, TOL);
rho_minus = evolve_rho_t(L_minus, rho_0, t, dim, TOL);
for n = 1:Nt
    F = calc_fidelity(rho_minus(:,:,n), rho_plus(:,:,n));
    QFI_om(n) = 2*(1 - F) / (1 + sqrt(F)) / d_om^2;
end

% 2. QFI for lambda
fprintf('Calculating QFI for lambda...\n');
L_plus = build_Lindbladian(om_base, lam_base+d_lam, ga_base, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base, lam_base-d_lam, ga_base, sx, sz, sm, I, Is, dim);
rho_plus = evolve_rho_t(L_plus, rho_0, t, dim, TOL);
rho_minus = evolve_rho_t(L_minus, rho_0, t, dim, TOL);
for n = 1:Nt
    F = calc_fidelity(rho_minus(:,:,n), rho_plus(:,:,n));
    QFI_lam(n) = 2*(1 - F) / (1 + sqrt(F)) / d_lam^2;
end

% 3. QFI for gamma
fprintf('Calculating QFI for gamma...\n');
L_plus = build_Lindbladian(om_base, lam_base, ga_base+d_ga, sx, sz, sm, I, Is, dim);
L_minus = build_Lindbladian(om_base, lam_base, ga_base-d_ga, sx, sz, sm, I, Is, dim);
rho_plus = evolve_rho_t(L_plus, rho_0, t, dim, TOL);
rho_minus = evolve_rho_t(L_minus, rho_0, t, dim, TOL);
for n = 1:Nt
    F = calc_fidelity(rho_minus(:,:,n), rho_plus(:,:,n));
    QFI_ga(n) = 2*(1 - F) / (1 + sqrt(F)) / d_ga^2;
end

% 4. QFI for time (quantum clock)
fprintf('Calculating QFI for time t...\n');
L_base = build_Lindbladian(om_base, lam_base, ga_base, sx, sz, sm, I, Is, dim);
rho_t_all = evolve_rho_t(L_base, rho_0, t, dim, TOL);

for n = 1:Nt
    if n == 1
        rho_curr = squeeze(rho_t_all(:,:,n));
        rho_next = squeeze(rho_t_all(:,:,n+1));
        F = calc_fidelity(rho_curr, rho_next);
        QFI_t(n) = 8 * (1 - sqrt(F)) / dt^2;
    elseif n == Nt
        rho_prev = squeeze(rho_t_all(:,:,n-1));
        rho_curr = squeeze(rho_t_all(:,:,n));
        F = calc_fidelity(rho_prev, rho_curr);
        QFI_t(n) = 8 * (1 - sqrt(F)) / dt^2;
    else
        rho_prev = squeeze(rho_t_all(:,:,n-1));
        rho_next = squeeze(rho_t_all(:,:,n+1));
        F = calc_fidelity(rho_prev, rho_next);
        QFI_t(n) = 2 * (1 - F) / (1 + sqrt(F)) / dt^2;
    end
end

fprintf('Computation done. \n');

%% Helper functions

function [R_sort, L_sort, lambda_sort] = sortingEigenvalues(dim, TOL, L)
    % Sort eigenvalues and eigenmatrices of Liouvillian L for time evolution.
    dim2 = dim*dim;
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
            if abs(real(an)-real(bm))<TOL && abs(imag(an)-imag(bm))<TOL && count<=dim2
                ind_RL(count, :) = [n, m];
                count = count + 1;
            end
        end
    end
    
    lambda = eig_R(ind_RL(:,1));
    [~, ind] = sort(real(lambda), 'descend');
    lambda_sort = lambda(ind);
    
    R_sort = cell(1, dim2);
    L_sort = cell(1, dim2);
    for k = 1:dim2
        Rk = reshape(R(:,ind_RL(ind(k),1)), dim, dim);
        Lk = reshape(Lmat(:,ind_RL(ind(k),2)), dim, dim);
        Ck = trace(Lk * Rk);
        if abs(Ck) > 1e-12
            Rk = Rk / sqrt(Ck);
            Lk = Lk / sqrt(Ck);
        end
        R_sort{k} = Rk;
        L_sort{k} = Lk;
    end
end

function L = build_Lindbladian(om, lam, ga, sx, sz, sm, I, Is, dim)
    % Build Lindbladian superoperator for two spins.
    H1 = om*kron(sz, I); 
    H2 = om*kron(I, sz); 
    HI = lam*kron(sx, sx);
    Hs = H1 + H2 + HI;
    
    L_H = -1i*(kron(Is, Hs) - kron(Hs.', Is));
    
    sm1 = kron(sm, I); 
    sm2 = kron(I, sm);
    L1 = ga*(kron(conj(sm1),sm1) -0.5*kron(Is,sm1'*sm1) -0.5*kron(sm1.'*conj(sm1),Is));
    L2 = ga*(kron(conj(sm2),sm2) -0.5*kron(Is,sm2'*sm2) -0.5*kron(sm2.'*conj(sm2),Is));
    
    L = L_H + L1 + L2;
end

function F = calc_fidelity(rho, sigma)
    % Fidelity for 2x2 density matrices (generalized for any dimension).
    M = rho * sigma;
    lambda = eig(M);
    lambda = real(lambda);
    lambda(lambda < 0) = 0;
    F = abs(sum(sqrt(lambda)))^2;
end

function rho_t = evolve_rho_t(L, rho_0, t_vec, dim, TOL)
    % Evolve initial state under Liouvillian L for all times in t_vec.
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
            rho = rho + c_k(k) * exp(lambda_sort(k)*t_vec(n)) * R_sort{k};
        end
        rho = (rho + rho')/2;
        rho = rho / trace(rho);
        rho_t(:,:,n) = rho;
    end
end