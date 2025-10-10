clc; clear; close all;

% 读取 Excel 数据
filename = '3.xlsx';  % 你的 Excel 文件名
data = xlsread(filename);     % 读取数据
centers = [];  % 用于存储圆心坐标
radii = [];    % 用于存储半径

% 解析 Excel 数据
for i = 1:3:84  % 每3列为一组 (x, y, r)
    x_i = data(:, i);     % 读取 x 坐标
    y_i = data(:, i + 1); % 读取 y 坐标
    r_i = data(:, i + 2); % 读取半径

    centers = [centers; x_i, y_i]; % 累加圆心
    radii = [radii; r_i];          % 累加半径
end

% 画图
figure; hold on; axis equal;
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');

% 画黑色圆圈（传感器覆盖范围）
viscircles(centers, radii, 'Color', 'k', 'LineWidth', 1);

% 画蓝色虚线边界框
plot([-1, 1, 1, -1, -1], [-1, -1, 1, 1, -1], 'b--', 'LineWidth', 2);

% 标题
title('EMT-PD-T1', 'Interpreter', 'latex');
box on;  % 保持四个轴显示

% 设置坐标轴范围
xlim([-1.2, 1.2]);
ylim([-1.2, 1.2]);
grid on;

% 保存图片为 PNG
% axis tight;
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.01));
% saveas(gcf, 'sensor_coverage.png','tight');
% exportgraphics(gcf, 'sensor_coverage.png', 'BackgroundColor', 'white');
print('EMT.png', '-dpng', '-r300');
