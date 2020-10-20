%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time."
% By Masakazu Emoto and Takeki Sunakawa
% This code calculates the bias-correlation to solve the law of motion by XPA algorithm
% This algoithm is reffered by Sunakawa (2020 Computational Economics)
% Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [zeta, psix, muxz, mpvec] = bias(gds, pds)
    
    global gamma rho alpha delta la intx x com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin intZ gridZ dZ ddZ
    global quada quadx quadK quadZ
    
    % Container
    muxz = zeros(intx,intZ); psix = zeros(intx,intZ);
    mpvec = zeros(intx,intZ); zeta = zeros(intx,intZ);
    
    % Calculate the bias-correlation 
    for iz = 1:intZ
        for ix = 1:intx
            % Population conditioned on labor productivity
            muxz(ix, iz) = sum(gds(:, ix)*da);
            % The amount of capital capital conditioned on labor productivity at deterministic steady state
            Know = sum(gds(:,ix)'*grida*da)/muxz(ix, iz);
            
            % Liner interpolation for policy function approximation at aggregate capital
            % Grid search for Know
            ia = gridsearch(Know, grida); 
            weight = (grida(ia + 1) - Know)/da;
            
            % Policy function at the amount of capital indexed to productivity
            mpvec(ix,iz) = (1 - weight) * pds(ia + 1, ix) + weight * pds(ia, ix);
            
            % Bias correlation
            zeta(ix, iz) = -(mpvec(ix,iz) - (sum(pds(:,ix)'*gds(:,ix)*da)/muxz(ix, iz)));
        end
        
        % Calculate the fraction of e-indexed capital (Sunakawa(2020 Computational Economics) page 28)
        mnow = 0; % Aggregate  wealth at deterministic steady state
        for ix = 1:intx
            psix(ix, iz) = sum(gds(:,ix)'*grida*da)/muxz(ix, iz);
            mnow = mnow + (sum(gds(:,ix)'*grida*da));
        end
    end
    % The ratio of the capital conditioned on labor productivity to aggregate capital
    psix = psix/mnow;
end