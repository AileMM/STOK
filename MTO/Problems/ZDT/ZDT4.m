function [Obj, Con] = ZDT4(var, M, opt, g)
[ps, D] = size(var);  % ps: number of solutions
f1 = var(:,1);
g_val = 1 + 10*(D - 1) + sum(var(:,2:end).^2 - 10*cos(4*pi*var(:,2:end)), 2);
h = 1 - sqrt(f1 ./ g_val);
f2 = g_val .* h;

Obj = [f1, f2];      % N¡Á2
Con = zeros(ps, 1);  % ÎÞÔ¼Êø
end

