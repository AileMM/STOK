function [Obj, Con] = DTLZ1(var, M, opt)
% DTLZ1 function
%   - var: design variable vector
%   - M: number of objectives
%   - opt: shift vector (for problem-specific shifting)
%   - g: Objective value shift

[ps, D] = size(var);

g  = 100*(D-M+1+sum((var(:,M:end)-0.5).^2-cos(20.*pi.*(var(:,M:end)-0.5)),2));
Obj = 0.5*repmat(1+g,1,M).*fliplr(cumprod([ones(size(var,1),1),var(:,1:M-1)],2)).*[ones(size(var,1),1),1-var(:,M-1:-1:1)];

% % Adding objective value shift (if needed)
% Obj = Obj + g;

% Constrain (if applicable, for example)
Con = zeros(ps, 1);

end
