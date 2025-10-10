% function model = buildMTGP(dec, obj, num_obj)
%    % Extract decision variables and objectives from the population
%    % %vertcat
%    pop_dec = dec;
%    pop_obj = obj;
%    
%    % Number of population samples
%    [num_pop, ~] = size(pop_dec);
%    
%    % Prepare training data for each objective
%    x_train = cell(1, num_obj);
%    y_train = cell(1, num_obj);
%    for i = 1:num_obj
%        x_train{i} = pop_dec;
%        y_train{i} = pop_obj(:, i);
%    end
%    
%    % Randomly initialize the interaction matrix C
%    C = rand(num_obj, num_obj) * 2 - 1;
%    
%    % Transform the objective values F to get H using H = C * F
%    H = zeros(num_obj, num_pop);
%    for i = 1:num_pop
%        F = [y_train{1}(i); y_train{2}(i)]; % For a two-objective case
%        H(:, i) = C * F;
%    end
%    
%    % Initialize model cell array for each transformed objective H
%    model = cell(1, num_obj);
%    
%    % Fit Gaussian Process models for each row of H 这里有转置，两个目标会得到两个模型
%    for t = 1:num_obj
%        model{t} = fitrgp(x_train{t}, H(t, :)', 'KernelFunction', 'squaredexponential');
%    end
   
   
function model = buildMTGP(dec, obj, num_obj)
    % dec: 决策变量矩阵
    % obj: 目标值矩阵
    % num_obj: 目标的数量
    
    % 提取决策变量和目标值
    pop_dec = dec;
    pop_obj = obj;
    
    % 获取种群数量
    [num_pop, ~] = size(pop_dec);
    
    % 准备训练数据
    x_train = cell(1, num_obj);
    y_train = cell(1, num_obj);
    for i = 1:num_obj
        x_train{i} = pop_dec;
        y_train{i} = pop_obj(:, i);
    end
    
    % 随机初始化交互矩阵 C 这个C需要有范围的吧，这里后续得看看论文
    C = rand(num_obj, num_obj) * 2;
    
    % 使用 H = C * F 转换目标值
    H = zeros(num_obj, num_pop);
    for i = 1:num_pop
        F = pop_obj(i, :)'; % 当前个体的目标值 (列向量)
        H(:, i) = C * F; % 计算 H
    end
    
    % 初始化用于存储每个目标 H 的模型
    model = cell(1, num_obj);
    
    % 为每个目标 H 训练高斯过程模型
    for t = 1:num_obj
        model{t} = fitrgp(x_train{t}, H(t, :)', 'KernelFunction', 'squaredexponential');
    end
end

