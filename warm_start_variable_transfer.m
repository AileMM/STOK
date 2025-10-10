function [database, FEsUsed_add, success, source] = warm_start_variable_transfer(database, knowledge_base, fun, M)
    success = false;
    FEsUsed_add = 0;
    try
        target_pop = database(:, 1:end-M);
        target_obj = database(:, end-M+1:end);
        dimt = size(target_pop, 2);
        source_pops = {knowledge_base.solutions};
        mind = cellfun(@(x) size(x, 2), source_pops);
        dims = min(mind);
        proj_dim = min(dimt, dims);

        % Step 1: Find most similar source task 这里是随机投影+MMD
        most_similar_idx_history = [];
        [most_similar_idx, ~] = FindMostSimilarTask(target_pop, source_pops, proj_dim);
        most_similar_idx_history(end+1) = most_similar_idx;   % 记录这次的索引
        disp(most_similar_idx_history);
        transfer_source = knowledge_base(most_similar_idx);
        idx = transfer_source.key_var;  %var_d
        transfer_var = mean(transfer_source.solutions(:, idx)); % 这里是取了均值的

        % Step 2: 从最差的帕累托前沿开始，逐层选择个体，直到达到替换数量
        pop_size = size(target_pop, 1);
        [~, FrontNo] = NSGASort(target_obj);
        max_front = max(FrontNo);
        %% =======这里有一个消融变量，即选择变量替换的比列======
        num_replace = ceil(0.2 * pop_size);
        selected_ind = [];  % 初始化选择个体下标
        for front = max_front:-1:1
            inds = find(FrontNo == front);
            selected_ind = [selected_ind; inds(:)];
            if length(selected_ind) >= num_replace
                selected_ind = selected_ind(1:num_replace);
                break;
            end
        end
        if length(selected_ind) < num_replace
            remain = setdiff(1:pop_size, selected_ind);
            extra = remain(randperm(length(remain), num_replace - length(selected_ind)));
            selected_ind = [selected_ind; extra(:)];
        end
        
        var = transfer_var;
        selected_idx = idx;
        % Step 2.5: 去除超过目标维度的变量
        valid_mask = selected_idx <= dimt;
        selected_idx = selected_idx(valid_mask);
        var = var(valid_mask);
        % 如果源任务的变量数小于目标任务的维度
        if length(var) < dimt
            % 不需要做任何改变，直接用源任务的变量
            selected_idx = idx(1:length(var));  % 根据源任务变量的数量限制索引
        end
        if length(var) > dimt
            % 如果源任务变量数大于目标维度，随机选择 dimt 个索引
            selected_idx = idx(randperm(length(idx), dimt));
            var = transfer_var(selected_idx);  % 确保替换的变量维度一致
        end

        % Step 5: 替换目标种群中的个体
        for i = 1:length(selected_ind)
            ind = selected_ind(i);
            % 检查 selected_idx 维度和 var 维度是否匹配
            if length(selected_idx) ~= length(var)
                error('Selected index dimensions do not match transfer variable length.');
            end
            % 进行替换，只替换目标种群中指定维度的变量
            target_pop(ind, selected_idx) = var;
            target_obj(ind, :) = fun(target_pop(ind, :));  % 重新计算目标值
        end

        % Step 4: Update database with modified individuals
        database = [target_pop, target_obj];
        FEsUsed_add = length(selected_ind); 
        success = true;
        source = transfer_source;
    catch
        % 异常时输出参数设为默认值
        database = [target_pop, target_obj];
        FEsUsed_add = 0;
        success = false;
        source = [];  
        return;
    end
end
