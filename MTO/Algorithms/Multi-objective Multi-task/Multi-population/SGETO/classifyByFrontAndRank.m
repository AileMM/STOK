% function [convergence_indices, diversity_indices] = classifyByFrontAndRank(front_no, rank, diversity_threshold)
%     % 基于 FrontNo 和 Rank 对个体进行分类
%     %
%     % 参数:
%     %   front_no: 个体的 FrontNo 列表（数组）
%     %   rank: 每个 FrontNo 下个体的 rank 值（数组）
%     %   diversity_threshold: Rank 的多样性分界比例（默认前后 10%）
%     %
%     % 返回:
%     %   convergence_indices: 收敛性特征强的个体索引
%     %   diversity_indices: 多样性特征强的个体索引
% 
%     if nargin < 3
%         diversity_threshold = 0.1; % 默认前后 10% 作为多样性边界
%     end
% 
%     num_individuals = length(front_no);
%     convergence_indices = [];
%     diversity_indices = [];
% 
%     % 获取唯一的 FrontNo
%     unique_fronts = unique(front_no);
% 
%     for i = 1:length(unique_fronts)
%         front = unique_fronts(i);
% 
%         % 找到属于该 FrontNo 的个体索引
%         front_indices = find(front_no == front);
% 
%         % 获取这些索引对应的 rank
%         front_ranks = rank(front_indices);
% 
%         % 按 rank 对个体进行排序
%         [~, sorted_order] = sort(front_ranks);
%         sorted_indices = front_indices(sorted_order);
% 
%         % 确定多样性特征的边界数量
%         boundary_count = max(1, floor(length(sorted_indices) * diversity_threshold));
% 
%         % 分类
%         if front <= 2
%             % FrontNo <= 2 的个体归为收敛性强   这里应该是往后接吧
%             convergence_indices = [convergence_indices, sorted_indices];
%         else
%             % FrontNo > 2 的个体需要进一步划分
%             diversity_set  = [diversity_indices, ...
%                                  sorted_indices(1:boundary_count), ...
%                                  sorted_indices(end-boundary_count+1:end)];
%                              
%             middle_set = sorted_indices(boundary_count+1:end-boundary_count);
% 
%             % 去除中间部分与边界的重叠（可能因边界数量较多导致重复）
%             middle_set = setdiff(middle_set, diversity_set);
% 
%             % 更新分类
%             diversity_indices = [diversity_indices, diversity_set];
%             convergence_indices = [convergence_indices, middle_set];
%             
% %             if length(sorted_indices) > 2 * boundary_count
% %                 convergence_indices = [convergence_indices, ...
% %                                        sorted_indices(boundary_count+1:end-boundary_count)];
% %             end
%         end
%     end
%     % 去重确保两个集合互斥
%     diversity_indices = setdiff(diversity_indices, convergence_indices);
% end
function [convergence_indices, diversity_indices] = classifyByFrontAndRank(front_no, rank, diversity_threshold)
    if nargin < 3
        diversity_threshold = 0.1;
    end

    num_individuals = length(front_no);
    convergence_indices = [];
    diversity_indices = [];

    unique_fronts = unique(front_no);

    for i = 1:length(unique_fronts)
        front = unique_fronts(i);

        front_indices = find(front_no == front);
        front_ranks = rank(front_indices);

        [~, sorted_order] = sort(front_ranks);
        sorted_indices = front_indices(sorted_order);

        boundary_count = max(1, floor(length(sorted_indices) * diversity_threshold));

        if front <= 2
            convergence_indices = [convergence_indices, sorted_indices]; % 改为竖向拼接
        else
            diversity_set = [sorted_indices(1:boundary_count), ...
                             sorted_indices(end-boundary_count+1:end)];
            middle_set = sorted_indices(boundary_count+1:end-boundary_count);

            middle_set = setdiff(middle_set, diversity_set); % 去重
            diversity_indices = [diversity_indices, diversity_set];
            convergence_indices = [convergence_indices, middle_set];
        end
    end

    diversity_indices = setdiff(diversity_indices, convergence_indices);
end
