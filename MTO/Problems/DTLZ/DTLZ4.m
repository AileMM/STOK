function [Obj, Con] = DTLZ4(var, M, opt)
% DTLZ4 function
%   - var: design variable vector
%   - M: number of objectives
%   - opt: shift vector (for problem-specific shifting)
%   - g: Objective value shift

[ps, D] = size(var);

var(:,1:M-1) = var(:,1:M-1).^100;
g  = sum((var(:,M:end)-0.5).^2,2);
Obj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(var(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(var(:,M-1:-1:1)*pi/2)];
% Constrain (if applicable, for evarample)
Con = zeros(ps, 1);

end
