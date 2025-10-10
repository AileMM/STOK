function [nondominated_set, nondominated_obj] = checkNonDominated(pop1, pop2)
    % pop1: 种群 A 的个体（带 Obj 属性的结构体数组）
    % pop2: 种群 B 的个体（带 Obj 属性的结构体数组）
    % 返回非支配解的集合：决策变量集合和对应目标值集合

    % 合并 pop1 和 pop2 的目标值和决策向量
    a_Obj = vertcat(pop1.Obj);  % A 中所有个体的目标值矩阵
    Obj = vertcat(pop2.Obj);    % B 中所有个体的目标值矩阵
    Dec = vertcat(pop2.Dec);    % B 中所有个体的决策变量矩阵
    
    num_B = size(Obj, 1);       % B 中个体数
    num_A = size(a_Obj, 1);     % A 中个体数

    is_nondominated = true(num_B, 1); % 初始化为全非支配解

    % 遍历 B 中的每个个体
    for i = 1:num_B
        b_obj = Obj(i, :); % B 中第 i 个体的目标值

        % 遍历 A 中每个个体是否支配 B 中个体
        for j = 1:num_A
            a = a_Obj(j, :); % A 中第 j 个体的目标值

            % 如果 A 中某个个体支配 B 中当前个体
            if all(a <= b_obj) && any(a < b_obj)
                is_nondominated(i) = false;  % 标记为非非支配解
                break;  % 跳过剩余比较
            end
        end
    end

    % 筛选非支配解
    nondominated_set = Dec(is_nondominated, :);  % 决策向量
    nondominated_obj = Obj(is_nondominated, :);  % 对应目标值
end
