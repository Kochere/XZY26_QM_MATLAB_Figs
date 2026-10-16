%% 计算差值并储存到三维矩阵中
%----------------------------------
t1 = 10;
t2 = 1000;
p0span = -0.5-1e-8: 2e-4: 0.5+1e-8; % -0.5<=p0<=0.5
f0span = 0-1e-8: 2e-4: 0.5+1e-8; % f0 的范围
% p0span = -0.5-1e-8: 1e-5: 0.5+1e-8; % -0.5<=p0<=0.5
% f0span = 0-1e-8: 1e-2: 0.5+1e-8; % f0 的范围
%----------------------------------

ga = 0.01; %gamma
ga0 = 0.01; %gamma0
om = 1; %beta
be = 1; %omega

np0= length(p0span);
nf0= length(f0span);

%----------------------------------
dd1 = zeros(np0, nf0); 
dd101 = zeros(np0, nf0);
dd2 = zeros(np0, nf0); 
%----------------------------------



% 计算 t1 和 t2 的结果
for ip = 1:np0
    p0 = p0span(ip);
    f0_max = max(f0span);  % 动态计算 f0 的最大允许值
    dd1(ip, :) = compute_dd(p0, f0_max, nf0, t1, be, ga, ga0, om);
    dd2(ip, :) = compute_dd(p0, f0_max, nf0, t2, be, ga, ga0, om);
end

% 定义一个辅助函数来计算 DD_func
function dd = compute_dd(p0, f0_max, nf0, t, be, ga, ga0, om)
    f0vals = linspace(0, f0_max, nf0);
    dd = zeros(1, nf0);
    for if0 = 1:nf0
        f0 = f0vals(if0);
        dd(if0) = DD_func(p0, f0, t, be, ga, ga0, om);
        if p0^2 + f0^2 > 0.25+1e-8
            dd(if0) = NaN;
        end
    end
end


