function optimum = W3_Optimum(N,M)
    X = linspace(0,1,N)';
    X = [X,zeros(N,M-2)+0.5,zeros(N,1)];
    optimum = linear(X);
    optimum = repmat(2:2:2*M,size(optimum,1),1).*optimum;
end

function Output = linear(x)
    Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end