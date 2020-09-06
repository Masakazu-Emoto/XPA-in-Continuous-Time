function [zeta, psix, muxz, mpvec] = bias(gds, pds)
    
% This code is to calculate the bias-correlation
    global gamma rho alpha delta la intx x com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin intZ gridZ ddZ
    
    for iz = 1:intZ
        for ix = 1:intx
            muxz(ix, iz) = sum(gds(:, ix)*da);
            Know = sum(gds(:,ix)'*grida*da)/muxz(ix, iz);
            
            % Liner interpolation for policy function approximation at aggregate wealth
            ia = gridsearch(Know, grida);
            weight = (grida(ia + 1) - Know)/da;

            mpvec(ix,iz) = (1 - weight) * pds(ia + 1, ix) + weight * pds(ia, ix);
            zeta(ix, iz) = -(mpvec(ix,iz) - (sum(pds(:,ix)'*gds(:,ix)*da)/muxz(ix, iz)));
        end
        
        % Calculate the fraction of e-indexed capital (Sunakawa(2019) page 28)
        mnow = 0; % Aggregate  wealth at deterministic steady state
        for ix = 1:intx
            psix(ix, iz) = sum(gds(:,ix)'*grida*da)/muxz(ix, iz);
            mnow = mnow + (sum(gds(:,ix)'*grida*da));
        end
    end
    psix = psix/mnow;
end