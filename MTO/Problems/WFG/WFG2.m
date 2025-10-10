function [Obj, Con] = WFG2(var, M, opt, g)
    % WFG2 function
    %   - var: design variable vector
    %   - M: number of objectives
    %   - opt: shift vector
    %   - g: Objective value move

    [N, D] = size(var);
    K = M-1;
    L = D - K;
    S = 2 : 2 : 2*M;
    A = ones(1,M-1);

    % 归一化到 [0,1]
    z01 = var ./ repmat(2:2:size(var,2)*2,N,1);

    % 第一阶段转换
    t1 = zeros(N,K+L);
    t1(:,1:K)     = z01(:,1:K);
    t1(:,K+1:end) = s_linear(z01(:,K+1:end), 0.35);

    % 第二阶段转换
    num_pairs = floor(L/2);
    t2 = zeros(N,K+num_pairs);
    t2(:,1:K) = t1(:,1:K);
    t2(:,K+1:K+num_pairs) = (t1(:,K+1:2:K+2*num_pairs-1) + ...
                             t1(:,K+2:2:K+2*num_pairs) + ...
                             2*abs(t1(:,K+1:2:K+2*num_pairs-1) - ...
                                   t1(:,K+2:2:K+2*num_pairs))) / 3;

    % 第三阶段转换
    t3 = zeros(N,M);
    for i = 1:M-1
        t3(:,i) = r_sum(t2(:,(i-1)*K/(M-1)+1:i*K/(M-1)), ones(1, K/(M-1)));
    end
    t3(:,M) = r_sum(t2(:,K+1:K+num_pairs), ones(1, num_pairs));

    % 计算目标值
    x = zeros(N,M);
    for i = 1:M-1
        x(:,i) = max(t3(:,M), A(:,i)) .* (t3(:,i) - 0.5) + 0.5;
    end
    x(:,M) = t3(:,M);

    h      = convex(x);
    h(:,M) = disc(x);
    Obj = repmat(x(:,M), 1, M) + repmat(S, N, 1) .* h;

    % 无约束
    Con = zeros(N, 1);
end

function Output = s_linear(y, A)
    Output = abs(y - A) ./ abs(floor(A - y) + A);
end

function Output = r_sum(y, w)
    Output = sum(y .* repmat(w, size(y,1),1),2) ./ sum(w);
end

function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1), 1-cos(x(:,1:end-1)*pi/2)],2)) .* ...
             [ones(size(x,1),1), 1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = disc(x)
    Output = 1 - x(:,1) .* (cos(5*pi*x(:,1))).^2;
end
