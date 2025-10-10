function distance = computeWasserstein(population1, population2)
    % 计算两个种群之间的 Wasserstein 距离 (EMD)
    % 输入：
    %   population1: 上一代种群 (N1 x D) 矩阵，其中 N1 为种群个体数，D 为维度
    %   population2: 当前代种群 (N2 x D) 矩阵，其中 N2 为种群个体数，D 为维度
    % 输出：
    %   distance: Wasserstein 距离（标量）

    % 确认输入维度相同
    assert(size(population1, 2) == size(population2, 2), '维度不匹配');

    % 计算两个种群之间的距离矩阵 (N1 x N2)
    N1 = size(population1, 1);
    N2 = size(population2, 1);
    costMatrix = zeros(N1, N2);

    for i = 1:N1
        for j = 1:N2
            % 计算欧几里得距离
            costMatrix(i, j) = norm(population1(i, :) - population2(j, :));
        end
    end

    % 初始化种群权重（均匀分布）
    weight1 = ones(N1, 1) / N1;
    weight2 = ones(N2, 1) / N2;

    % 计算 Wasserstein 距离 (使用线性规划)
    distance = emd_mex(weight1, weight2, costMatrix);
end

function distance = emd_mex(p1, p2, costMatrix)
    % 使用 MATLAB 的线性规划函数来最小化运输成本
    f = costMatrix(:);
    Aeq = [kron(eye(size(p1, 1)), ones(1, size(p2, 1))); ...
           kron(ones(1, size(p1, 1)), eye(size(p2, 1)))];
    beq = [p1; p2];
    
    % 线性规划求解最小化成本
    options = optimoptions('linprog', 'Display', 'none');
    flow = linprog(f, [], [], Aeq, beq, zeros(size(f)), [], options);
    
    % 计算总成本即为 Wasserstein 距离
    distance = sum(f .* flow);
end
