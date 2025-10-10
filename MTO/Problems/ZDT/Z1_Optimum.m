function optimum = Z1_Optimum(N,M)
    optimum(:,1) = linspace(0,1,N)';
    optimum(:,2) = 1 - optimum(:,1).^0.5;
end