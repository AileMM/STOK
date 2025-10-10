function optimum = W4_Optimum(N,M)
    optimum = UniformPoint(N,M);
    optimum = optimum./repmat(sqrt(sum(optimum.^2,2)),1,M);
    optimum = repmat(2:2:2*M,size(optimum,1),1).*optimum;
end