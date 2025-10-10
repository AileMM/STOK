%% Step 1: 读取数据
% 假设数据在 Excel 的第1个sheet里，从第1行开始
data = xlsread('迁移变量变化.xlsx', 1);   % data 是 10×4261 矩阵

% 假设下标在 Excel 的第2个sheet里，从第1列开始
idx = xlsread('迁移变量变化.xlsx', 2);    % idx 是 32×1 向量

%% Step 2: 提取对应下标的数据
selected_data = data(:, idx);

%% Step 3: 转换为 table 并输出为新的 Excel
T = array2table(selected_data);

% 给列名，直接用 idx 值
varNames = "Var_" + string(idx(:)');
T.Properties.VariableNames = varNames;

% 写入新的 Excel 文件
writetable(T, 'after_selected_data.xlsx');

disp('处理完成！已生成 selected_data.xlsx，列名已对应下标');
