function optimum = Z3_Optimum(N,M)
    optimum(:,1) = linspace(0,1,N)';
    optimum(:,2) = 1 - optimum(:,1).^0.5 - optimum(:,1).*sin(10*pi*optimum(:,1));
    optimum      = optimum(NDSort(optimum,1)==1,:);
end