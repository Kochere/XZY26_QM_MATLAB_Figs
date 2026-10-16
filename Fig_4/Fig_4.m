%% Fig_4.m
% Two-qubit dissipative model: compare two observables for estimating lambda.
% Computes Var(lambda), Our bound, and QCRB for two observables.
% Saves workspace variables: t, Var1, Bound1, QCRB1, Var2, Bound2, QCRB2.

clear; close all;

%% Physical parameters and grid
dim = 4;
Is = eye(dim); 
I  = eye(2);

sx = [0 1;1 0];  sy = [0 -1i;1i 0];  sz = [1 0;0 -1];
sm = [0 0;1 0];

om_base = 0.5;   lam_base = 1.0;   ga_base = 0.2;
TOL = 1e-10;
d_om = 0.001;   d_lam = 0.001;   d_ga = 0.001;

% Initial state |↓↓>
down = [0; 1];
psi_0 = kron(down, down);
rho_0 = psi_0 * psi_0';

% Time grid
t = linspace(1e-2, 30+1e-2, 6000);
Nt = length(t);
dt = t(2) - t(1);

% Two observables
Op1 = om_base * kron(sz, I);   % observable 1
Op2 = kron(sx, sx);            % observable 2

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
    for ii = 1:dim2
        an = eig_R(ii);
        for jj = 1:dim2
            bm = eig_L(jj);
            if abs(real(an)-real(bm))<TOL && abs(imag(an)-imag(bm))<TOL && count<=dim2
                ind_RL(count, :) = [ii, jj];
                count = count + 1;
            end
        end
    end
    lambda = eig_R(ind_RL(:,1));
    [~, ind] = sort(real(lambda), 'descend');
    lambda_sort = lambda(ind);
    R_sort = cell(1, dim2);
    L_sort = cell(1, dim2);
    for ii = 1:dim2
        Rk = reshape(R(:, ind_RL(ind(ii),1)), dim, dim);
        Lk = reshape(Lmat(:, ind_RL(ind(ii),2)), dim, dim);
        Ck = trace(Lk * Rk);
        if abs(Ck) > 1e-12
            Rk = Rk / sqrt(Ck);
            Lk = Lk / sqrt(Ck);
        end
        R_sort{ii} = Rk;
        L_sort{ii} = Lk;
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
    for ii = 1:Nc
        c_k(ii) = trace(rho_0 * L_sort{ii});
    end
    Nt = length(t_vec);
    rho_t = zeros(dim, dim, Nt);
    for ii = 1:Nt
        rho = zeros(dim, dim);
        for kk = 1:Nc
            rho = rho + c_k(kk) * exp(lambda_sort(kk) * t_vec(ii)) * R_sort{kk};
        end
        rho = (rho + rho') / 2;
        rho = rho / trace(rho);
        rho_t(:, :, ii) = rho;
    end
end

%% Generic function to compute for a given observable
function [Var, Bound, QCRB] = compute_for_obs(Op, rho_0, om_base, lam_base, ga_base, ...
        d_om, d_lam, d_ga, t, dt, Nt, dim, TOL, sx, sz, sm, I, Is)
    
    % ---------- QFI (independent of observable) ----------
    % omega
    L_p = build_Lindbladian(om_base+d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
    L_m = build_Lindbladian(om_base-d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
    rho_p = evolve_rho(L_p, rho_0, t, dim, TOL);
    rho_m = evolve_rho(L_m, rho_0, t, dim, TOL);
    QFI_om = zeros(1,Nt);
    for ii = 1:Nt
        F = calc_fidelity(rho_m(:,:,ii), rho_p(:,:,ii));
        QFI_om(ii) = 2*(1-sqrt(F))/d_om^2;
    end
    QFI_om = max(QFI_om, 0);
    
    % lambda
    L_p = build_Lindbladian(om_base, lam_base+d_lam, ga_base, sx, sz, sm, I, Is, dim);
    L_m = build_Lindbladian(om_base, lam_base-d_lam, ga_base, sx, sz, sm, I, Is, dim);
    rho_p = evolve_rho(L_p, rho_0, t, dim, TOL);
    rho_m = evolve_rho(L_m, rho_0, t, dim, TOL);
    QFI_lam = zeros(1,Nt);
    for ii = 1:Nt
        F = calc_fidelity(rho_m(:,:,ii), rho_p(:,:,ii));
        QFI_lam(ii) = 2*(1-sqrt(F))/d_lam^2;
    end
    QFI_lam = max(QFI_lam, 0);
    QCRB = 1 ./ QFI_lam;
    
    % gamma
    L_p = build_Lindbladian(om_base, lam_base, ga_base+d_ga, sx, sz, sm, I, Is, dim);
    L_m = build_Lindbladian(om_base, lam_base, ga_base-d_ga, sx, sz, sm, I, Is, dim);
    rho_p = evolve_rho(L_p, rho_0, t, dim, TOL);
    rho_m = evolve_rho(L_m, rho_0, t, dim, TOL);
    QFI_ga = zeros(1,Nt);
    for ii = 1:Nt
        F = calc_fidelity(rho_m(:,:,ii), rho_p(:,:,ii));
        QFI_ga(ii) = 2*(1-sqrt(F))/d_ga^2;
    end
    QFI_ga = max(QFI_ga, 0);
    
    % time
    L_base = build_Lindbladian(om_base, lam_base, ga_base, sx, sz, sm, I, Is, dim);
    rho_base = evolve_rho(L_base, rho_0, t, dim, TOL);
    QFI_t = zeros(1,Nt);
    for ii = 1:Nt
        if ii == 1
            F = calc_fidelity(rho_base(:,:,ii), rho_base(:,:,ii+1));
            QFI_t(ii) = 8*(1-sqrt(F))/dt^2;
        elseif ii == Nt
            F = calc_fidelity(rho_base(:,:,ii-1), rho_base(:,:,ii));
            QFI_t(ii) = 8*(1-sqrt(F))/dt^2;
        else
            F = calc_fidelity(rho_base(:,:,ii-1), rho_base(:,:,ii+1));
            QFI_t(ii) = 2*(1-F)/(1+sqrt(F))/dt^2;
        end
    end
    QFI_t = max(QFI_t, 0);
    
    % ---------- Observable expectation and derivatives ----------
    exp_H = zeros(1,Nt); exp_H2 = zeros(1,Nt);
    for ii = 1:Nt
        exp_H(ii) = real(trace(Op * rho_base(:,:,ii)));
        exp_H2(ii) = real(trace(Op * Op * rho_base(:,:,ii)));
    end
    Var_O = exp_H2 - exp_H.^2;
    
    % d/dλ
    L_p = build_Lindbladian(om_base, lam_base+d_lam, ga_base, sx, sz, sm, I, Is, dim);
    L_m = build_Lindbladian(om_base, lam_base-d_lam, ga_base, sx, sz, sm, I, Is, dim);
    rho_p = evolve_rho(L_p, rho_0, t, dim, TOL);
    rho_m = evolve_rho(L_m, rho_0, t, dim, TOL);
    exp_p = zeros(1,Nt); exp_m = zeros(1,Nt);
    for ii = 1:Nt
        exp_p(ii) = real(trace(Op * rho_p(:,:,ii)));
        exp_m(ii) = real(trace(Op * rho_m(:,:,ii)));
    end
    dH_dLam = (exp_p - exp_m) / (2*d_lam);
    
    % d/dω
    L_p = build_Lindbladian(om_base+d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
    L_m = build_Lindbladian(om_base-d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
    rho_p = evolve_rho(L_p, rho_0, t, dim, TOL);
    rho_m = evolve_rho(L_m, rho_0, t, dim, TOL);
    for ii = 1:Nt
        exp_p(ii) = real(trace(Op * rho_p(:,:,ii)));
        exp_m(ii) = real(trace(Op * rho_m(:,:,ii)));
    end
    dH_dOm = (exp_p - exp_m) / (2*d_om);
    
    % d/dγ
    L_p = build_Lindbladian(om_base, lam_base, ga_base+d_ga, sx, sz, sm, I, Is, dim);
    L_m = build_Lindbladian(om_base, lam_base, ga_base-d_ga, sx, sz, sm, I, Is, dim);
    rho_p = evolve_rho(L_p, rho_0, t, dim, TOL);
    rho_m = evolve_rho(L_m, rho_0, t, dim, TOL);
    for ii = 1:Nt
        exp_p(ii) = real(trace(Op * rho_p(:,:,ii)));
        exp_m(ii) = real(trace(Op * rho_m(:,:,ii)));
    end
    dH_dGa = (exp_p - exp_m) / (2*d_ga);
    
    % d/dt
    dH_dt = zeros(1,Nt);
    dH_dt(1) = (exp_H(2) - exp_H(1))/dt;
    dH_dt(end) = (exp_H(end) - exp_H(end-1))/dt;
    dH_dt(2:end-1) = (exp_H(3:end) - exp_H(1:end-2))/(2*dt);
    
    % ---------- Var and Our bound ----------
    Var = Var_O ./ (dH_dLam.^2 + 1e-12);
    
    denom = abs(dH_dLam);
    denom(denom < 1e-15) = 1e-15;
    ratio_lam = ones(size(t));
    ratio_om  = abs(dH_dOm ./ denom).^2;
    ratio_ga  = abs(dH_dGa ./ denom).^2;
    ratio_t   = abs(dH_dt ./ denom).^2;
    
    minQ = 1e-15;
    QFI_lam_safe = QFI_lam; QFI_lam_safe(QFI_lam_safe < minQ) = inf;
    QFI_om_safe  = QFI_om;  QFI_om_safe(QFI_om_safe < minQ)   = inf;
    QFI_ga_safe  = QFI_ga;  QFI_ga_safe(QFI_ga_safe < minQ)   = inf;
    QFI_t_safe   = QFI_t;   QFI_t_safe(QFI_t_safe < minQ)     = inf;
    
    quantity_lam = ratio_lam ./ QFI_lam_safe;
    quantity_om  = ratio_om  ./ QFI_om_safe;
    quantity_ga  = ratio_ga  ./ QFI_ga_safe;
    quantity_t   = ratio_t   ./ QFI_t_safe;
    
    Bound = max([quantity_lam; quantity_om; quantity_ga; quantity_t], [], 1);
    Bound(~isfinite(Bound)) = 1e-20;
    
    maxY = 1e12;
    Var(Var > maxY)   = maxY;
    Bound(Bound > maxY) = maxY;
    QCRB(QCRB > maxY) = maxY;
end

%% Compute for both observables
fprintf('Computing for observable 1: omega*sz⊗I ...\n');
[Var1, Bound1, QCRB1] = compute_for_obs(Op1, rho_0, om_base, lam_base, ga_base, ...
    d_om, d_lam, d_ga, t, dt, Nt, dim, TOL, sx, sz, sm, I, Is);

fprintf('Computing for observable 2: sx⊗sx ...\n');
[Var2, Bound2, QCRB2] = compute_for_obs(Op2, rho_0, om_base, lam_base, ga_base, ...
    d_om, d_lam, d_ga, t, dt, Nt, dim, TOL, sx, sz, sm, I, Is);

fprintf('Done. Variables ready: t, Var1, Bound1, QCRB1, Var2, Bound2, QCRB2.\n');