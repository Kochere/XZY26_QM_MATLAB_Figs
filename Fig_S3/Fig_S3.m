%% Fig_S3.m
% Compute derivatives of <H^(1)> of t, time, omega, lambda, gamma,
% and verify the scaling law for an open two-spin system.
% Saves variables: t, dO_dt, dO_dOm, dO_dLam, dO_dGa, LHS, RHS.

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

% Finite-difference steps
d_om = 0.001;
d_lam = 0.001; 
d_ga = 0.001;

% Initial state |↓↓>
down = [0; 1];
psi_0 = kron(down, down);
rho_0 = psi_0 * psi_0';

% Time grid
t = linspace(0.01, 30, 90000);
Nt = length(t);
dt = t(2)-t(1);

% Observable: H^(1) = om*(sz⊗I + I⊗sz) + lam*(sx⊗sx)
Op = om_base*kron(sz, I) + lam_base*kron(sx, sx);  % Note: om_base and lam_base used

%% Base evolution
L_base = build_L(om_base, lam_base, ga_base, sx, sz, sm, I, Is, dim);
[vals_base, ~, purity_base, conc_base] = evolve_all(L_base, rho_0, t, dim, TOL, Op);
O_vals = vals_base;   % expectation of Op

%% Compute derivatives w.r.t. parameters using central differences
% d/d_om
L_plus = build_L(om_base+d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
L_minus = build_L(om_base-d_om, lam_base, ga_base, sx, sz, sm, I, Is, dim);
[vals_p, ~, ~, ~] = evolve_all(L_plus, rho_0, t, dim, TOL, Op);
[vals_m, ~, ~, ~] = evolve_all(L_minus, rho_0, t, dim, TOL, Op);
dO_dOm = (vals_p - vals_m)/(2*d_om);

% d/d_lam
L_plus = build_L(om_base, lam_base+d_lam, ga_base, sx, sz, sm, I, Is, dim);
L_minus = build_L(om_base, lam_base-d_lam, ga_base, sx, sz, sm, I, Is, dim);
[vals_p, ~, ~, ~] = evolve_all(L_plus, rho_0, t, dim, TOL, Op);
[vals_m, ~, ~, ~] = evolve_all(L_minus, rho_0, t, dim, TOL, Op);
dO_dLam = (vals_p - vals_m)/(2*d_lam);

% d/d_ga
L_plus = build_L(om_base, lam_base, ga_base+d_ga, sx, sz, sm, I, Is, dim);
L_minus = build_L(om_base, lam_base, ga_base-d_ga, sx, sz, sm, I, Is, dim);
[vals_p, ~, ~, ~] = evolve_all(L_plus, rho_0, t, dim, TOL, Op);
[vals_m, ~, ~, ~] = evolve_all(L_minus, rho_0, t, dim, TOL, Op);
dO_dGa = (vals_p - vals_m)/(2*d_ga);

% d/dt (numerical)
dO_dt = zeros(size(O_vals));
dO_dt(1) = (O_vals(2) - O_vals(1))/dt;
dO_dt(end) = (O_vals(end) - O_vals(end-1))/dt;
for n = 2:Nt-1
    dO_dt(n) = (O_vals(n+1) - O_vals(n-1))/(2*dt);
end

% Scaling law quantities
LHS = t .* dO_dt;
RHS = om_base*dO_dOm + lam_base*dO_dLam + ga_base*dO_dGa;

fprintf('Computation done. \n');

%% Helper functions (unchanged)
function [R_sort, L_sort, lambda_sort] = sortingEigenvalues(dim, TOL, L)
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

function L = build_L(om, lam, ga, sx, sz, sm, I, Is, dim)
    H1 = om*kron(sz,I); 
    H2 = om*kron(I,sz); 
    HI = lam*kron(sx,sx);
    Hs = H1+H2+HI;
    L_H = -1i*(kron(Is, Hs) - kron(Hs.', Is));
    
    sm1 = kron(sm, I); 
    sm2 = kron(I, sm);
    L1 = ga*(kron(conj(sm1),sm1) -0.5*kron(Is,sm1'*sm1) -0.5*kron(sm1.'*conj(sm1),Is));
    L2 = ga*(kron(conj(sm2),sm2) -0.5*kron(Is,sm2'*sm2) -0.5*kron(sm2.'*conj(sm2),Is));
    L = L_H + L1 + L2;
end

function [vals, rho_all, purity_vals, concurrence_vals] = evolve_all(L, rho_0, t_vec, dim, TOL, Op)
    [R_sort, L_sort, lambda_sort] = sortingEigenvalues(dim, TOL, L);
    Nc = length(lambda_sort);
    c_k = zeros(1, Nc);
    for k = 1:Nc
        c_k(k) = trace(rho_0 * L_sort{k});
    end
    
    Nt = length(t_vec);
    vals = zeros(1, Nt);
    rho_all = zeros(dim, dim, Nt);
    purity_vals = zeros(1, Nt);
    concurrence_vals = zeros(1, Nt);
    
    for n = 1:Nt
        rho = zeros(dim, dim);
        for k = 1:Nc
            rho = rho + c_k(k) * exp(lambda_sort(k)*t_vec(n)) * R_sort{k};
        end
        rho = (rho + rho')/2;
        rho = rho / trace(rho);
        rho_all(:,:,n) = rho;
        
        if ~isempty(Op)
            vals(n) = real(trace(Op * rho));
        end
        
        purity_vals(n) = real(trace(rho*rho));
        
        sy = [0 -1i; 1i 0];
        R = rho * kron(sy, sy) * conj(rho) * kron(sy, sy);
        lambda_eig = sort(real(eig(R)), 'descend');
        concurrence_vals(n) = max(0, sqrt(lambda_eig(1)) - sqrt(lambda_eig(2)) - sqrt(lambda_eig(3)) - sqrt(lambda_eig(4)));
    end
end