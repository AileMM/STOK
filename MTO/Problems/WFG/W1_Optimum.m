function optimum = W1_Optimum(N,M)
   
    optimum = UniformPoint(N,M);
    c = ones(size(optimum,1),M);
    for i = 1 : size(optimum,1) 
        for j = 2 : M
            temp = optimum(i,j)/optimum(i,1)*prod(1-c(i,M-j+2:M-1));
            c(i,M-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
        end
    end
    x = acos(c)*2/pi;
    temp = (1-sin(pi/2*x(:,2))).*optimum(:,M)./optimum(:,M-1);
    a = 0 : 0.0001 : 1;
    E = abs(temp*(1-cos(pi/2*a))-1+repmat(a+cos(10*pi*a+pi/2)/10/pi,size(x,1),1));
    [~,rank] = sort(E,2);
    for i = 1 : size(x,1)
        x(i,1) = a(min(rank(i,1:10)));
    end
    optimum      = convex(x);
    optimum(:,M) = mixed(x);
    optimum      = repmat(2:2:2*M,size(optimum,1),1).*optimum;
end

function Output = convex(x)
    Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = mixed(x)
    Output = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
end