function [Obj, Con] = WFG3(var, M, opt, g)
% WFG3 function
%   - var: design variable vector
%   - M: number of objectives
%   - opt: shift vector
%   - g: Objective value move
[N,D] = size(var);

K = M-1;
L = D - K;  
D = 1;
S = 2 : 2 : 2*M;
A = [1,zeros(1,M-2)];

% z01 = var./repmat(2:2:size(var,2)*2,N,1);
% t1 = zeros(N,K+L);
% t1(:,1:K)     = z01(:,1:K);
% t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);
% t2 = zeros(N,floor(K+L/2));
% t2(:,1:K) = t1(:,1:K);
% % Same as <t2(:,i)=r_nonsep(t1(:,K+2*(i-K)-1:K+2*(i-K)),2)>
% t2(:,K+1:floor(K+L/2)) = (t1(:,K+1:2:end) + t1(:,K+2:2:end) + 2*abs(t1(:,K+1:2:end)-t1(:,K+2:2:end)))/3;
% % ---------------------------------------------------------
% t3 = zeros(N,M);
% for i = 1 : M-1
%     t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
% end
% t3(:,M) = r_sum(t2(:,K+1:K+L/2),ones(1,L/2));
% 
% x = zeros(N,M);
% for i = 1 : M-1
%     x(:,i) = max(t3(:,M),A(:,i)).*(t3(:,i)-0.5)+0.5;
% end

% 第一步：t1 构造
Z = var ./ repmat(2:2:(size(var,2)*2), N, 1);  % 生成除数，逐元素除法
t1 = zeros(N, K + L);
t1(:, 1:K) = Z(:, 1:K);
t1(:, K+1:end) = s_linear(Z(:, K+1:end), 0.35);  % 自定义函数
% 第二步：t2 构造
L_half = floor(L / 2);
t2 = zeros(N, K + L_half);
t2(:, 1:K) = t1(:, 1:K);
% 奇偶对组合：(K+1)::2 和 (K+2)::2
idx1 = K + (1:2:L);
idx2 = K + (2:2:L);
len_pairs = min(length(idx1), length(idx2));
for i = 1:len_pairs
    t2(:, K+i) = (t1(:, idx1(i)) + t1(:, idx2(i)) + 2 * abs(t1(:, idx1(i)) - t1(:, idx2(i)))) / 3;
end
% 第三步：t3 构造
t3 = ones(N, M);
% 计算每组的大小
K_divide_M_sub_1 = floor(K / (M - 1));
% 前 M-1 个维度，每个维度聚合 K_divide_M_sub_1 个特征
for i = 1:(M - 1)
    range_idx = ((i - 1) * K_divide_M_sub_1 + 1):(i * K_divide_M_sub_1);
    t3(:, i) = r_sum(t2(:, range_idx), ones(1, K_divide_M_sub_1));
end

% 最后一维使用 t2 的后 L/2 个元素
t3(:, M) = r_sum(t2(:, K+1:end), ones(1, L_half));
x(:,M) = t3(:,M);
% ========================================================

h   = linear(x);
Obj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

% No constraint in WFG3
Con = zeros(N, 1);
end

function Output = s_linear(y,A)
    Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = r_nonsep(y,A)
    Output = zeros(size(y,1),1);
    for j = 1 : size(y,2)
        Temp = zeros(size(y,1),1);
        for k = 0 : A-2
            Temp = Temp+abs(y(:,j)-y(:,1+mod(j+k,size(y,2))));
        end
        Output = Output+y(:,j)+Temp;
    end
    Output = Output./(size(y,2)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = linear(x)
    Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end
