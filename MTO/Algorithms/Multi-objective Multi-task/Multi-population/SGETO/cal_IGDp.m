function result = cal_IGDp(MTOData, varargin)
% <Metric> <Multi-objective>
row = 1;
results = struct(); 
resultsDec = struct(); 

for t = 1:MTOData.Problems.T
    fieldName = sprintf('Task%d', t); % 动态生成字段名
    results.(fieldName) = []; % 初始化为空数组
    resultsDec.(fieldName) = [];
end

    for prob = 1:length(MTOData.Problems)
        optimum = MTOData.Problems(prob).Optimum;
        for task = 1:MTOData.Problems(prob).T
            for algo = 1:length(MTOData.Algorithms)
                gen = size(MTOData.Results(prob, algo, 1).Obj{task}, 1);
                igdp = zeros(MTOData.Reps, gen);
                BestObj = {};
                %如果没有搞错的话 这里应该只有任务1的IGD值，改过来了
                for rep = 1:MTOData.Reps
                    for g = 1:gen
                        Obj = squeeze(MTOData.Results(prob, algo, rep).Obj{task}(g, :, :));
                        Dec = squeeze(MTOData.Results(prob, algo, rep).Dec{task}(g, :, :));
                        CV = squeeze(MTOData.Results(prob, algo, rep).CV(task, g, :));
%                         BestObj{rep} = getBestObj(Obj, CV);
                       for i = 1:length(Obj)
                           obj = Obj{i}; %获取每一代的目标函数值
                           dec = Dec{i}; %获取每一代的决策变量
                           outdec{rep, i} = dec;
                           igdp(rep, i) = getIGDp(obj, optimum{task}); 
                      
                       end
                        
                    end
                end
%                 if task == 1
%                     igdp_task1 = igdp;  % 存储任务1的结果
%                     results(task).Task1 = igdp_task1;
%                 elseif task == 2
%                     igdp_task2 = igdp;  % 存储任务2的结果
%                     results(task).Task2 = igdp_task2;
%                 end
              
                %% 这里保留两个任务的优化结果的IGD值 最后输出即可
%                 result.TableData(row, algo, :) = igdp(:, end);
                % 指定文件名

%                 for rep = 1:MTOData.Reps
%                     result.ConvergeData.Y(row, algo, rep, :) = igdp(rep, :);
%                     result.ConvergeData.X(row, algo, rep, :) = [1:gen] ./ gen .* MTOData.Problems(prob).maxFE ./ MTOData.Problems(prob).T;
%                     result.ParetoData.Obj{row, algo, rep} = squeeze(BestObj{rep}(:, :));
%                 end
            end
            result.ParetoData.Optimum{row}(:, :) = optimum{task};
            row = row + 1;
%             for task = 1:MTOData.Problems.T
            fieldName = sprintf('Task%d', t);
            results(task).(fieldName) = igdp; % 动态存储任务结果
     
            resultsDec(task).(fieldName) = outdec; 
%             end
        end
         % 生成保存的文件名（包含算法名称）
        output_filename = sprintf('%s_%s_IGDp10.mat', MTOData.Algorithms.Name, MTOData.Problems(1).Name);
        output_dec = sprintf('%s_%s_Dec10.mat', MTOData.Algorithms.Name, MTOData.Problems(1).Name);

        % 保存结果 评价指标和决策变量
        save(output_filename, 'results');
        save(output_dec, 'resultsDec');
     end
end
 
function BestObj = getBestObj(Obj, CV)
Feasible = find(all(CV <= 0, 2));
if isempty(Feasible)
    Best = [];
else
    Best = NDSort(Obj(Feasible, :), 1) == 1;
end
BestObj = Obj(Feasible(Best), :);
end
 

