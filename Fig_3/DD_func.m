function value = DD_func(p0, f0, t, be, ga, ga0, om)
 
    % 计算热浴参数
    N = 1.0 / (exp(be * om) - 1.0);  % Bose-Einstein 分布
    gamma_up = ga * N;
    gamma_down = ga * (N + 1);
    
    gamma_p = gamma_up + gamma_down;
    gamma_m = gamma_up - gamma_down;
    
    % 稳态值
    p_steady = gamma_m / (2 * gamma_p);
    
    % 时间演化
    p_t = p_steady + (p0 - p_steady) * exp(-gamma_p * t);
    f_t = f0 * exp(-1j * om * t) * exp(-(2 * ga0 + gamma_p / 2) * t);
    f_t_mag_sq = abs(f_t)^2;
    
    %本征值
    muup = 1./2 + sqrt(p_t.^2 + f_t_mag_sq);
    mudown = 1./2 - sqrt(p_t.^2 + f_t_mag_sq);

    H =  (exp(- gamma_p * t) + gamma_p * t -1)/ (gamma_p*t);
    T = 1/be;

    % term of Quantum Fisher Informtion of beta
    A = om^2 * t^2 * T^4 *gamma_up^2 *gamma_down^2 /(gamma_p^2 * p_steady^2 *muup*mudown);
    B = p_steady * (1-H) * f_t_mag_sq /(p_t - p_steady);
    C = 4*p_t*(p_t - p_steady*H);
    D = (1/4 -p_t^2)* (2-(p_steady*(H-1))/(p_t-p_steady));


    value = A * B * (C+D);
end

