function [Obj, Con] = ZDT3(var, M, opt, g)
[ps, D] = size(var);  % ps: number of solutions
f1 = var(:,1);  % 第一目标
g_val = 1 + 9 * mean(var(:,2:end), 2);  % g(x)
h = 1 - sqrt(f1 ./ g_val) - (f1 ./ g_val) .* sin(10 * pi * f1);  % h(x)
f2 = g_val .* h;  % 第二目标
Obj = [f1, f2];  % N×2

Con = zeros(ps, 1);  % 无约束
end

