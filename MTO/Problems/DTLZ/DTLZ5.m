function [Obj, Con] = DTLZ5(var, M, opt)
% DTLZ5 function
%   - var: design variable vector
%   - M: number of objectives
%   - opt: shift vector (for problem-specific shifting)
%   - g: Objective value shift

[ps, D] = size(var);

g    = sum((var(:,M:end)-0.5).^2,2);
Temp = repmat(g,1,M-2);
var(:,2:M-1) = (1+2*Temp.*var(:,2:M-1))./(2+2*Temp);
Obj = repmat(1+g,1,M).*fliplr(cumprod([ones(size(g,1),1),cos(var(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(var(:,M-1:-1:1)*pi/2)];

Con = zeros(ps, 1);

end
