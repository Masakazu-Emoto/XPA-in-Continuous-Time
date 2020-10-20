%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time."
% By Masakazu Emoto and Takeki Sunakawa
% This code calculates new forecastig rule by XPA algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Kdotnew, xpavec] = outer(ps, muxz, psix, zeta)

    global gamma rho alpha delta la intx x com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin intZ gridZ dZ ddZ
    global quada quadx quadK quadZ
    
    % Container
    xpavec = zeros(intx,intZ); Kdotnew = zeros(intK,intZ); bias = 1;
    
    % Calculate New law of motion
    for ik = 1:intK
        for iz = 1:intZ
            
            for ix = 1:intx
                % The amount of capital conditioned on labor productivity
                Know = psix(ix, iz)* gridK(ik);
                
                % Grid search for Know
                ia = gridsearch(Know,grida);
                weight = (grida(ia + 1) - Know)/da; 
                
                % Liner interpolation for policy function approximation at aggregate capital
                gvec = (1 - weight) * ps(ia + 1, ix, ik, iz) + weight * ps(ia, ix, ik, iz);
                if bias == 1
                    xpavec(ix,iz) = gvec + zeta(ix,iz); % with bias-correlation
                else
                    xpavec(ix,iz) = gvec; % without bias-correlation
                end
            end
            
            % New law of motion
            Kdotnew(ik, iz) =  sum(sum(muxz(:,iz)' * xpavec(:,iz),1));
        end
    end
end