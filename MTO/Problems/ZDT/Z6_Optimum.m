function optimum = Z6_Optimum(N,M)
    minf1  = 0.280775;
    optimum(:,1) = linspace(minf1,1,N)';
    optimum(:,2) = 1 - optimum(:,1).^2;
end