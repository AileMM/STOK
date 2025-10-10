function [Obj, Con] = DTLZ7(var, M, opt, g)
% DTLZ7 function
%   - var: design variable vector
%   - M: number of objectives
%   - opt: shift vector (for problem-specific shifting)
%   - g: Objective value shift

[ps, D] = size(var);
Obj = zeros(size(var,1),M);
g = 1+9*mean(var(:,M:end),2);
Obj(:,1:M-1) = var(:,1:M-1);
Obj(:,M) = (1+g).*(M-sum(Obj(:,1:M-1)./(1+repmat(g,1,M-1)).*(1+sin(3*pi.*Obj(:,1:M-1))),2));

% Constraints (if applicable)
Con = zeros(ps, 1);

end
