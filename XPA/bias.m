function [muxz, psix, zeta] = bias(gds, pds)
%% bias.m : This code calculates the bias-correlation terms to solve for the forecasting rule by XPA algorithm
%
%% INPUTS
%    gds        : the stationary distribution
%    pds        : the policy function for savings at the deterministic steady state
%% OUTPUTS
%
%    muxz       : the population conditioned on labor productivity
%    psix       : the ratio of the capital conditioned on labor productivity to aggregate capital
%    zeta       : bias correction terms

    global gamma rho alpha delta la intx x com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin intZ gridZ dZ ddZ
    global quada quadx quadK quadZ
    
    % Container
    muxz = zeros(intx,intZ); psix = zeros(intx,intZ);
    mpvec = zeros(intx,intZ); zeta = zeros(intx,intZ);
    
    % Calculate the bias-correlation terms
    for iz = 1:intZ
        for ix = 1:intx
            % the population conditioned on labor productivity
            muxz(ix, iz) = sum(gds(:, ix)*da);
            % the amount of capital conditioned on labor productivity at the deterministic steady state
            Know = sum(gds(:,ix)'*grida*da)/muxz(ix, iz);
            
            % Liner interpolation for policy function approximation at the conditional capital
            % Grid search for Know
            ia = gridsearch(Know, grida); 
            weight = (grida(ia + 1) - Know)/da;
            mpvec(ix,iz) = (1 - weight) * pds(ia + 1, ix) + weight * pds(ia, ix);
            
            % Bias correction
            zeta(ix, iz) = -(mpvec(ix,iz) - (sum(pds(:,ix)'*gds(:,ix)*da)/muxz(ix, iz)));
        end
        
        % Calculate the fraction of e-indexed capital (Sunakawa (2020 Computational Economics) page 28)
        mnow = 0; % Aggregate wealth at the deterministic steady state
        for ix = 1:intx
            psix(ix, iz) = sum(gds(:,ix)'*grida*da)/muxz(ix, iz);
            mnow = mnow + (sum(gds(:,ix)'*grida*da));
        end
    end
    % The ratio of the capital conditioned on labor productivity to aggregate capital
    psix = psix/mnow;
    
end