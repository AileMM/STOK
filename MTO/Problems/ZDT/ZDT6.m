function [Obj, Con] = ZDT6(var, M, opt, g)
[ps, D] = size(var);
f1 = 1 - exp(-4 * var(:,1)) .* (sin(6 * pi * var(:,1))).^6;
g_val = 1 + 9 * (mean(var(:,2:end), 2)).^0.25;
h = 1 - (f1 ./ g_val).^2;
f2 = g_val .* h;

Obj = [f1, f2];      % N¡Á2
Con = zeros(ps, 1);  % ÎÞÔ¼Êø
end
