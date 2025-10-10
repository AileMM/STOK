% % 读取Excel数据
% data = readtable('plot bar.xlsx');
% 
% % 提取变量名和数值
% varNames    = data.kcat;
% beforeVals  = data.Before_transfer;
% afterVals   = data.After_transfer;
% 
% % 构造堆叠数据
% diffVals = afterVals - beforeVals;   % 增量（可能正也可能负）
% Y = [beforeVals, diffVals];          % 两部分堆叠
% 
% % 绘制堆叠柱状图
% figure;
% b = bar(Y,'stacked','BarWidth',0.6);
% 
% % 设置颜色
% b(1).FaceColor = [0.3 0.6 0.9];  % 蓝色 (Before)
% b(2).FaceColor = [0.9 0.4 0.4];  % 红色 (After-Before)
% 
% % 坐标轴设置
% set(gca, 'XTick', 1:length(varNames), ...
%          'XTickLabel', varNames, ...
%          'FontName','Times New Roman', ...
%          'FontSize', 12);
% xtickangle(45); % 变量名倾斜，防止重叠
% 
% xlabel('k_{cat} Variables', 'FontSize',14,'FontName','Times New Roman');
% ylabel('k_{cat} Value', 'FontSize',14,'FontName','Times New Roman');
% title('k_{cat} Before and After Transfer (Stacked)', ...
%       'FontSize',16,'FontName','Times New Roman');
% 
% legend({'Before transfer','After-Before'}, ...
%        'FontSize',12,'FontName','Times New Roman','Location','best');
% 
% grid on;



% 读取Excel数据
data = readtable('plot bar.xlsx');

% 提取变量名和数值
varNames    = data.kcat;
beforeVals  = data.Before_transfer;
afterVals   = data.After_transfer;

% 构造堆叠数据
diffVals = afterVals - beforeVals;   % 增量
Y = [beforeVals, diffVals];          % 两部分堆叠

figure;

%% --- 上半部分：堆叠柱状图（log scale） ---
subplot(2,1,1);
b = bar(Y,'stacked','BarWidth',0.6);

% 设置颜色
% b(1).FaceColor = [0.4 0.76 0.92];  % 浅蓝 (Before)
% b(2).FaceColor = [1.0 0.60 0.40];  % 橙色 (After-Before)
b(1).FaceColor = [0.20 0.60 0.95];  % 更亮的蓝 (Before)
b(2).FaceColor = [1.00 0.40 0.20];  % 更鲜的橙 (After-Before)

% y轴对数坐标
set(gca,'YScale','log');

% 坐标轴设置

set(gca, 'XTick', 1:length(varNames), ...
         'XTickLabel', varNames, ...
         'FontName','Times New Roman', ...
         'FontSize', 10);
% xtickangle(45);
xticks(1:4:size(post,2));

ylabel('$k_{cat}$ Value (log scale)','Interpreter','latex','FontName','Times New Roman','FontSize',12);
title('$k_{cat}$ Before and After Transfer (Stacked)', ...
      'FontSize',14,'FontName','Times New Roman','Interpreter','latex');

legend({'Before transfer','After-Before'}, ...
       'FontSize',10,'FontName','Times New Roman','Location','NorthWest');

grid on;

%% --- 下半部分：变化比值 After / Before ---
subplot(2,1,2);
changeRatio = afterVals ./ beforeVals;  % 比值

bar(changeRatio,'FaceColor',[0.8 0.6 0.8]); % 马卡龙紫
% bar(changeRatio,'FaceColor',[0.85 0.35 0.95]); % 高饱和度紫
set(gca, 'XTick', 1:length(varNames), ...
         'XTickLabel', varNames, ...
         'FontName','Times New Roman', ...
         'FontSize', 10);
xticks(1:4:size(post,2));
% xtickangle(45);

ylabel('After / Before', 'FontSize',12,'FontName','Times New Roman');
title('Fold Change of $k_{cat}$ After vs. Before', ...
      'Interpreter','latex','FontSize',14,'FontName','Times New Roman');

grid on;

% 在柱子上标注数值
for i = 1:length(changeRatio)
    text(i, changeRatio(i), sprintf('%.2f', changeRatio(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',8, 'FontName','Times New Roman');
end
print(gcf,'plotbar.png','-dpng','-r600');

