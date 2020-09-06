function [Kdotnew, xpavec] = outer(ps, muxz, psix, zeta)
    
    global gamma rho alpha delta la intx x com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin intZ gridZ dZ ddZ
    global quada quadx quadK quadZ
    
    % Calculate New Forecasting Rules
    xpavec = zeros(intx,intZ); Kdotnew = zeros(intK,intZ);
    bias = 1;
    for ik = 1:intK
        for iz = 1:intZ
            
            for ix = 1:intx
                Know = psix(ix, iz)* gridK(ik);
                ia = gridsearch(Know,grida);
                weight = (grida(ia + 1) - Know)/da; 
                
                gvec = (1 - weight) * ps(ia + 1, ix, ik, iz) + weight * ps(ia, ix, ik, iz);
                if bias == 1
                    xpavec(ix,iz) = gvec + zeta(ix,iz); % with bias-correlation
                else
                    xpavec(ix,iz) = gvec; % without bias-correlation
                end
            end
            
            Kdotnew(ik, iz) =  sum(sum(muxz(:,iz)' * xpavec(:,iz),1));
        end
    end
end