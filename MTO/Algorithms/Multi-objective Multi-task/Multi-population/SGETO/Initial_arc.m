function Archive = Initial_arc(Algo, Prob, Individual_Class, varargin)
%% Multi-task - Initialize and evaluate the population
% Input: Algorithm, Problem, Individual_Class
% Output: population

%------------------------------- Copyright --------------------------------
% Copyright (c) Yanchi Li. You are free to use the MToP for research
% purposes. All publications which use this platform should acknowledge
% the use of "MToP" or "MTO-Platform" and cite as "Y. Li, W. Gong, F. Ming,
% T. Zhang, S. Li, and Q. Gu, MToP: A MATLAB Optimization Platform for
% Evolutionary Multitasking, 2023, arXiv:2312.08134"
%--------------------------------------------------------------------------

n = numel(varargin);
if n == 0
    N = Prob.N;
elseif n == 1
    N = varargin{1};
else
    return;
end

for t = 1:Prob.T
    for i = 1:N
        Archive{t}(i) = Individual_Class();
        lhs_sample = lhsdesign(1, max(Prob.D));
%         lhs_sample = lhsdesign(1, Prob.D(t)); %实际应用应该根据每个任务的维度来
%         Archive{t}(i).Dec = Prob.Lb{t} + lhs_sample .* (Prob.Ub{t} - Prob.Lb{t});
        Archive{t}(i).Dec = Prob.Lb{Prob.T} + lhs_sample .* (Prob.Ub{Prob.T} - Prob.Lb{Prob.T});

    end
    Archive{t} = Algo.Evaluation(Archive{t}, Prob, t);
end
end