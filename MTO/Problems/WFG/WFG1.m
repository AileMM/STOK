function [Obj, Con] = WFG1(var, M, opt, g)
% WFG1 function
%   - var: design variable vector
%   - M: number of objectives
%   - opt: shift vector
%   - g: Objective value move

[ps,D] = size(var);
K = M-1;
L = D - K;
D = 1;
S = 2 : 2 : 2*M;
A = ones(1,M-1);

z01 = var./repmat(2:2:size(var,2)*2,ps,1);

t1 = zeros(ps,K+L);
t1(:,1:K)     = z01(:,1:K);
t1(:,K+1:end) = s_linear(z01(:,K+1:end),0.35);

t2 = zeros(ps,K+L);
t2(:,1:K)     = t1(:,1:K);
t2(:,K+1:end) = b_flat(t1(:,K+1:end),0.8,0.75,0.85);

t3 = zeros(ps,K+L);
t3 = b_poly(t2,0.02);

t4 = zeros(ps,M);
for i = 1 : M-1
    t4(:,i) = r_sum(t3(:,(i-1)*K/(M-1)+1:i*K/(M-1)),2*((i-1)*K/(M-1)+1):2:2*i*K/(M-1));
end
t4(:,M) = r_sum(t3(:,K+1:K+L),2*(K+1):2:2*(K+L));

x = zeros(ps,M);
for i = 1 : M-1
    x(:,i) = max(t4(:,M),A(i)).*(t4(:,i)-0.5)+0.5;
end
x(:,M) = t4(:,M);

h      = convex(x);
h(:,M) = mixed(x);
Obj = repmat(D*x(:,M),1,M) + repmat(S,ps,1).*h;

% pso constraint in WFG1
Con = zeros(ps, 1);
end

function Output = s_linear(y,A)
    Output = abs(y-A)./abs(floor(A-y)+A);
end

function Output = b_flat(y,A,B,C)
    Output = A+min(0,floor(y-B))*A.*(B-y)/B-min(0,floor(C-y))*(1-A).*(y-C)/(1-C);
    Output = round(Output*1e4)/1e4;
end

function Output = b_poly(y,a)
    Output = y.^a;
end

function Output = r_sum(y,w)
    Output = sum(y.*repmat(w,size(y,1),1),2)./sum(w);
end

function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
    Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end
      

