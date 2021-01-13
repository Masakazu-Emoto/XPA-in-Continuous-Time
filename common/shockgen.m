function [Zsim,Zup,Zdown,zweight] = shockgen(N,mu,sigma,dT,dZ,Zmean,Zmin,Zmax,gridZ,shock)

    mmu = -1 + mu; 

    Zsim = zeros(N,1);
    for time = 1:N-1
        if time == 1
            Zsim(time+1) = mu * dT * Zmean + (1 - mu * dT) * Zmean + sigma * shock(time) * sqrt(dT);
%             Zsim(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT));
        else
            Zsim(time+1) = mu * dT * Zmean + (1 - mu * dT) * Zsim(time) + sigma * shock(time) * sqrt(dT);
%             Zsim(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT));
        end
    end

    Zup     = zeros(N,1);
    Zdown   = zeros(N,1);
    zweight = zeros(N,1);

    for time = 1:N
        Zsim(time)    = max([Zsim(time) Zmin+0.000001]);
        Zsim(time)    = min([Zsim(time) Zmax-0.000001]);
        Zdown(time)   = floor((Zsim(time) - Zmin)/dZ) + 1;
        Zup(time)     = ceil((Zsim(time) - Zmin)/dZ) + 1;
        zweight(time) = (gridZ(Zup(time)) - Zsim(time))/dZ;
    end

end