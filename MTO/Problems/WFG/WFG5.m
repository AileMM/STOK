function [Obj, Con] = WFG5(var, M, opt, g)
% WFG5 function
%   - var: design variable vector
%   - M: number of objectives
%   - opt: shift vector
%   - g: Objective value move

[N,D] = size(var);

K = M-1;
L = D - K;
D = 1;
S = 2 : 2 : 2*M;
A = ones(1,M-1);

z01 = var./repmat(2:2:size(var,2)*2,N,1);

t1 = zeros(N,K+L);
t1 = s_decept(z01,0.35,0.001,0.05);

t2 = zeros(N,M);
for i = 1 : M-1
    t2(:,i) = r_sum(t1(:,(i-1)*K/(M-1)+1:i*K/(M-1)),ones(1,K/(M-1)));
end
t2(:,M) = r_sum(t1(:,K+1:K+L),ones(1,L));

x = zeros(N,M);
for i = 1 : M-1
    x(:,i) = max(t2(:,M),A(:,i)).*(t2(:,i)-0.5)+0.5;
end
x(:,M) = t2(:,M);
h = concave(x);
Obj = repmat(D*x(:,M),1,M) + repmat(S,N,1).*h;

% No constraint in WFG5
Con = zeros(N, 1);
end

function Output = s_decept(y,A,B,C)
    Output = 1+(abs(y-A)-B).*(floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B);
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = concave(x)
    Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
