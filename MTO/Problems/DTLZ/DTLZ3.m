function [Obj, Con] = DTLZ3(var, M, opt)
% DTLZ3 function
%   - var: design variable vector
%   - M: number of objectives
%   - opt: shift vector (for problem-specific shifting)
%   - g: Objective value shift

[ps, D] = size(var);

% Calculate g
h = 100 * (D - M + 1 + sum((var(:, M:end) - 0.5).^2 - cos(20 * pi * (var(:, M:end) - 0.5)), 2));
% Objective calculation
Obj = repmat(1 + h, 1, M) .* fliplr(cumprod([ones(ps, 1), cos(var(:, 1:M-1) * pi / 2)], 2)) .* [ones(ps, 1), sin(var(:, M-1:-1:1) * pi / 2)];

% Constrain (if applicable, for example)
Con = zeros(ps, 1);

end
