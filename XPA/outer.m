function Kdotnew = outer(ps, muxz, psix, zeta)
%% outer.m : This code calculates new forecasting rule by XPA algorithm
%
%% INPUTS
%    ps         : the policy function for savings
%    muxz       : the population conditioned on labor productivity
%    psix       : the ratio of the capital conditioned on labor productivity to aggregate capital
%    zeta       : bias correction terms
%
%% OUTPUTS
%    Kdotnew    : new forecasting rule

    global gamma rho alpha delta la intx x com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin intZ gridZ dZ ddZ
    global quada quadx quadK quadZ
    
    % Container
    xpavec = zeros(intx,intZ); Kdotnew = zeros(intK,intZ); 
    % flag for bias correction
    bias = 1;
    
    % Calculate the new forecasting rule
    for ik = 1:intK
        for iz = 1:intZ
            
            for ix = 1:intx
                % the amount of capital conditioned on labor productivity
                Know = psix(ix, iz)* gridK(ik);
                
                % Liner interpolation for policy function approximation at the conditional capital
                % Grid search for Know
                ia = gridsearch(Know,grida);
                weight = (grida(ia + 1) - Know)/da;                
                gvec = (1 - weight) * ps(ia + 1, ix, ik, iz) + weight * ps(ia, ix, ik, iz);
                if bias == 1
                    xpavec(ix,iz) = gvec + zeta(ix,iz); % with bias-correlation
                else
                    xpavec(ix,iz) = gvec; % without bias-correlation
                end
            end
            
            % New forecasting rule
            Kdotnew(ik, iz) =  sum(sum(muxz(:,iz)' * xpavec(:,iz),1));
        end
    end
    
end