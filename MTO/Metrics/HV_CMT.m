function result = HV_CMT(MTOData, varargin)
% <Metric> <Multi-objective>

% Competitive Hypervolume of All Tasks

%------------------------------- Reference --------------------------------
% @Article{Li2022CompetitiveMTO,
%   author     = {Li, Genghui and Zhang, Qingfu and Wang, Zhenkun},
%   journal    = {IEEE Transactions on Evolutionary Computation},
%   title      = {Evolutionary Competitive Multitasking Optimization},
%   year       = {2022},
%   pages      = {1-1},
%   doi        = {10.1109/TEVC.2022.3141819},
% }
%--------------------------------------------------------------------------

%------------------------------- Copyright --------------------------------
% Copyright (c) Yanchi Li. You are free to use the MToP for research
% purposes. All publications which use this platform should acknowledge
% the use of "MToP" or "MTO-Platform" and cite as "Y. Li, W. Gong, F. Ming,
% T. Zhang, S. Li, and Q. Gu, MToP: A MATLAB Optimization Platform for
% Evolutionary Multitasking, 2023, arXiv:2312.08134"
%--------------------------------------------------------------------------

Par_flag = false;
if length(varargin) >= 1
    % Parallel Calculate
    Par_flag = varargin{1};
end
%============================
results = struct(); 
for t = 1:MTOData.Problems.T
    fieldName = sprintf('Task%d', t); % 动态生成字段名
    results.(fieldName) = []; % 初始化为空数组
end
%============================
result.Metric = 'Max';
result.RowName = {};
result.ColumnName = {};
% Data for Table
result.TableData = [];
% Data for Converge Plot
result.ConvergeData.X = [];
result.ConvergeData.Y = [];
% Data for Pareto Plot
result.ParetoData.Obj = {};
result.ParetoData.Optimum = [];

for prob = 1:length(MTOData.Problems)
    if MTOData.Problems(prob).M <= 1
        return;
    end
end
result.RowName = {MTOData.Problems.Name};
result.ColumnName = {MTOData.Algorithms.Name};

% Calculate Competitive HV
for prob = 1:length(MTOData.Problems)
    % Get Optimum
    AllBestObj = [];
    AllBestCV = [];
    for algo = 1:length(MTOData.Algorithms)
        for task = 1:MTOData.Problems(prob).T
            for rep = 1:MTOData.Reps
                AllBestObj = [AllBestObj; squeeze(MTOData.Results(prob, algo, rep).Obj{task}(end, :, :))];
%                 AllBestCV = [AllBestCV; squeeze(MTOData.Results(prob, algo, rep).CV(task, end, :))];
                AllBestCV = [AllBestCV; squeeze(MTOData.Results(prob, algo, rep).CV{task}(end, :, :))];
            end
        end
    end
%     optimum = getBestObj(AllBestObj, AllBestCV); %源代码所有结果的最优值作为参考点
    optimum = MTOData.Problems(prob).Optimum;
    for algo = 1:length(MTOData.Algorithms)
%         gen = size(MTOData.Results(prob, algo, 1).Obj{1}, 1); %源代码
        gen = size(MTOData.Results(prob, algo, 1).Obj{task}, 1);
        hv = zeros(MTOData.Reps, gen);
        BestObj = {};
        if Par_flag
            parfor rep = 1:MTOData.Reps
                for g = 1:gen
                    Obj = []; CV = [];
                    for task = 1:MTOData.Problems(prob).T
                        Obj_t = squeeze(MTOData.Results(prob, algo, rep).Obj{task}(g, :, :));
%                         CV_t = squeeze(MTOData.Results(prob, algo,
%                         rep).CV(task, g, :)); %源代码
                        CV_t = squeeze(MTOData.Results(prob, algo, rep).CV{task}(g, :, :));
                        Obj = [Obj; Obj_t];
                        CV = [CV; CV_t];
                    end
                    BestObj{rep} = getBestObj(Obj, CV);
                    hv(rep, g) = getHV(BestObj{rep}, optimum);
                    hv(rep, g) = getHV(BestObj{rep}, MTOData.Problems.Optimum);
%                     %========自己代码保存结果==============================
%                     for t = 1:size(Obj, 1)
%                         for f = 1:size(Obj, 2)
%                             if ~isempty(Obj{t, f})
%                                 hv(t, f) = getHV(Obj{t, f},MTOData.Problems.Optimum{t});  % 假设 calcHV 接收目标值
%                             else
%                                 hv(t, f) = 0;  % 没有可行解时HV为0
%                             end
%                         end
%                         %====================================
%                         fieldName = sprintf('Task%d', t);
%                         results(t).(fieldName) = hv(t, :); % 动态存储任务结果 
%                     end
%                     %======================================================
%                     output_filename = sprintf('%s_%s_HV.mat', MTOData.Algorithms.Name, MTOData.Problems(1).Name);
%                     % 保存结果 评价指标和决策变量
%                     save(output_filename, 'results');
                end
            end
        else
            for rep = 1:MTOData.Reps
                for g = 1:gen
                    Obj = []; CV = [];
                    for task = 1:MTOData.Problems(prob).T
                        Obj_t = squeeze(MTOData.Results(prob, algo, rep).Obj{task}(g, :, :));
%                         CV_t = squeeze(MTOData.Results(prob, algo, rep).CV(task, g, :));
                        CV_t = squeeze(MTOData.Results(prob, algo, rep).CV{task}(g, :, :));
                        Obj = [Obj; Obj_t];
                        CV = [CV; CV_t];
                        for i = 1:length(Obj_t)
                           obj = Obj_t{i}; %获取每一代的目标函数值
                           hv(rep, i) = getHV(obj, optimum{task});      
                        end
                       fieldName = sprintf('Task%d', task);
                       results(task).(fieldName) = hv; % 动态存储任务结果
                    end
                    BestObj{rep} = getBestObj(Obj, CV);
                end
            end
            
        end
         % 生成保存的文件名（包含算法名称）
        output_filename = sprintf('%s_%s_HV.mat', MTOData.Algorithms.Name, MTOData.Problems(1).Name);

        % 保存结果 评价指标和决策变量
        save(output_filename, 'results');
        end
        result.TableData(prob, algo, :) = hv(:, end);
        for rep = 1:MTOData.Reps
            result.ConvergeData.Y(prob, algo, rep, :) = hv(rep, :);
            result.ConvergeData.X(prob, algo, rep, :) = [1:gen] ./ gen .* MTOData.Problems(prob).maxFE;
            result.ParetoData.Obj{prob, algo, rep} = squeeze(BestObj{rep}(:, :));
        end
        result.ParetoData.Optimum{prob}(:, :) = optimum;
        
    end
end

% % 源码
% function BestObj = getBestObj(Obj, CV)
% Feasible = find(all(CV <= 0, 2));
% if isempty(Feasible)
%     Best = [];
% else
%     Best = NDSort(Obj(Feasible, :), 1) == 1;
% end
% BestObj = Obj(Feasible(Best), :);
% end

%===================自己代码====================
function BestObj = getBestObj(Obj, CV)
% Obj 和 CV 是 3x7 的 cell，每个 cell 是个体数 x M 的矩阵

[numTasks, numGens] = size(Obj);
BestObj = cell(numTasks, numGens);  % 初始化输出

for t = 1:numTasks
    for g = 1:numGens
        objMat = Obj{t, g};  % 当前任务第g代的目标函数值
        cvMat = CV{t, g};    % 当前任务第g代的约束值

        % 找到可行解（所有约束值 <= 0）
        feasibleIdx = all(cvMat <= 0, 2);

        if any(feasibleIdx)
            feasibleObj = objMat(feasibleIdx, :);
            rank = NDSort(feasibleObj, 1);  % 非支配排序
            bestIdx = (rank == 1);
            BestObj{t, g} = feasibleObj(bestIdx, :);  % 记录第一前沿解
        else
            BestObj{t, g} = [];  % 没有可行解
        end
    end
end

end




