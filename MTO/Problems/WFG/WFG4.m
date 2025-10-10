function [Obj, Con] = WFG4(var, M, opt, g)
% WFG4 function
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
t1 = s_multi(z01,30,10,0.35);

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

% No constraint in WFG4
Con = zeros(N, 1);
end

function Output = s_multi(y,A,B,C)
    Output = (1+cos((4*A+2)*pi*(0.5-abs(y-C)/2./(floor(C-y)+C)))+4*B*(abs(y-C)/2./(floor(C-y)+C)).^2)/(B+2);
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = concave(x)
    Output = fliplr(cumprod([ones(size(x,1),1),sin(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),cos(x(:,end-1:-1:1)*pi/2)];
end
