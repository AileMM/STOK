% function [Obj, Con] = ZDT1(var, M, opt, g)
% % WFG1 function
% %   - var: design variable vector
% %   - M: number of objectives
% %   - opt: shift vector
% %   - g: Objective value move
% [ps,D] = size(var);
% g = 1 + 9 * sum(var(2:end)) / (D - 1);
% f(1) = var(1);
% f(2) = g * (1 - sqrt(var(1) / g));
% Obj = f;
% 
% % pso constraint in WFG1
% Con = zeros(ps, 1);
% end

function [Obj, Con] = ZDT1(var, M, opt, g)

[ps, D] = size(var);  % ps: number of samples
g = 1 + 9 * sum(var(:, 2:end), 2) / (D - 1);  % N¡Á1
f1 = var(:, 1);                               % N¡Á1
f2 = g .* (1 - sqrt(f1 ./ g));                % N¡Á1

Obj = [f1, f2];                               % N¡Á2
Con = zeros(ps, 1);                           % N¡Á1, no constraint
end


