function [Obj, Con] = ZDT2(var, M, opt, g)
[ps, D] = size(var);  % ps: population size (number of individuals)
% Objective function calculations
f1 = var(:,1);  % First objective
g_val = 1 + 9 * mean(var(:,2:end), 2);  % g(x)
h = 1 - (f1 ./ g_val).^2;  % h(x)
f2 = g_val .* h;  % Second objective

% Output
Obj = [f1, f2];  % N ¡Á 2
Con = zeros(ps, 1);  % No constraints
end


