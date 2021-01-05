function [rds, wds, Kds, Ads, uds, cds, pds, ids, Vds, gds, X, Y, Z] = steadystate()
%% steadystate.m : This code calculates the deterministic steady state
% Reference : Achdou et al. (2017 NBER)
%% INPUTS : NA
%% OUTPUTS : 

%% Summary of the algorithm
% Step 0 : Initial guesses
% Step 1 : Solves for the deterministic steady state
% Step 1-1 : Value function iteration by the finite differential method
% Step 1-2 : Calculates the stationary distribution
% Step 1-3 : Checks the stationary equilibrium condition
% Step 2 : Calculates values at the steady state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global gamma rho alpha delta la intx x com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    
    %--------------------------------------------------%
    % Step 0 : Initial guesses
    %--------------------------------------------------%
    % Initial guess for r, K and w 
    r = 0.005; rmax = rho; rmin = 0.0001;
    K = (((alpha) / (r + delta)) ^ (1 / (1 - alpha))) * LAve;
    w = (1 - alpha) * (K ^ alpha) * ((LAve) ^ (-alpha));
    
    % Finite difference approximation of the partial derivatives
    Vaf = zeros(inta,intx); Vab = zeros(inta,intx);

    % Consumption
    c = zeros(inta,intx);
    
    % Initial guess for the value function
    if gamma == 1
        v0(:,1) = log((w * (1 - tau).* x(1) + w * com.* (1 - x(1)) + r.*grida))/rho;
        v0(:,2) = log((w * (1 - tau).* x(2) + w * com.* (1 - x(2)) + r.*grida))/rho;
    else
        v0(:,1) = (w * (1 - tau).* x(1) + w * com.* (1 - x(1)) + r.*grida).^(1-gamma)/(1-gamma)/rho;
        v0(:,2) = (w * (1 - tau).* x(2) + w * com.* (1 - x(2)) + r.*grida).^(1-gamma)/(1-gamma)/rho;
    end
    v = v0;

    %--------------------------------------------------%
    % Step 1 : Solves for the deterministic steady state
    %--------------------------------------------------%
    for iter=1:maxitK % Outer loop
%     disp('Main loop iteration')
%     disp(iter)
    
        %--------------------------------------------------%
        % Step 1-1 : Value function iteration by the finite differencial method
        %--------------------------------------------------%
        for n = 1:maxit % Inner loop
        
            V = v;

            % Forward Difference
            Vaf(1:inta-1,:) = (V(2:inta,:)-V(1:inta-1,:))/da;
            Vaf(inta,:) = (w * (1 - tau).* x + w * com.* (1 - x) + r.*amax).^(-gamma); %will never be used, but impose state constraint a<=amax just in case
   
            % Backward Difference
            Vab(2:inta,:) = (V(2:inta,:)-V(1:inta-1,:))/da;
            Vab(1,:) = (w * (1 - tau).* x + w * com.* (1 - x) + r.*amin).^(-gamma);  %state constraint boundary condition
    
            % Consumption and savings with forward difference
            cf = (Vaf).^(-1/gamma);
            sf = w * (1 - tau).* xx + w * com.* (1 - xx) + r.*aa - cf;
            
            % Consumption and savings with backward difference
            cb = (Vab).^(-1/gamma);
            sb = w * (1 - tau).* xx + w * com.* (1 - xx) + r.*aa - cb;
            
            % Consumption and derivative of value function at steady state
            c0 = w * (1 - tau).* xx + w * com.* (1 - xx)  + r.*aa;
            Va0 = c0.^(-gamma);

            % dV_upwind makes a choice of forward or backward differences based on
            % the sign of the drift
            If = sf > 0;              % Positive drift --> Forward difference
            Ib = sb < 0;              % Negative drift --> Backward difference
            I0 = (1 - If - Ib);       % Drift is zero --> At steady state

            Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; % important to include third term

            c = (Va_Upwind).^(-1/gamma);
            if gamma == 1
                u = log(c);
            else
                u = (c.^(1-gamma) )/(1-gamma);
            end
            
            % Construct sparse matrix for individual wealth
            X = -min(sb,0)/da;
            Y = -max(sf,0)/da + min(sb,0)/da;
            Z = max(sf,0)/da;

            A1 = spdiags(Y(:,1),0,inta,inta) + spdiags(X(2:inta,1),-1,inta,inta) + spdiags([0;Z(1:inta-1,1)],1,inta,inta);
            A2 = spdiags(Y(:,2),0,inta,inta) + spdiags(X(2:inta,2),-1,inta,inta) + spdiags([0;Z(1:inta-1,2)],1,inta,inta);
            Ads = [A1,sparse(inta,inta);sparse(inta,inta),A2] + Aswitch;
            
            B = (1/Delta + rho)*speye(inta*intx) - Ads;

            u_stacked = reshape(u,inta*intx,1);
            V_stacked = reshape(V,inta*intx,1);
            
            b = u_stacked + V_stacked/Delta;

            V_stacked = B\b; % SOLVE SYSTEM OF EQUATIONS
            V = reshape(V_stacked,inta,intx);

            Vchange = V - v;
            v = V;

            dist(n) = max(max(abs(Vchange)));
%             disp(max(max(abs(Vchange))))
            if dist(n) < crit
%                 disp('Value Function Converged, Iteration = ')
%                 disp(n)
                break
            end
            
        end % Inner loop

        %--------------------------------------------------%
        % Step 1-2 : Calculates the stationary distribution
        %--------------------------------------------------%
        AT = Ads';
        b = zeros(inta*intx,1);

        % need to fix one value, otherwise matrix is singular
        i_fix = 1;
        b(i_fix)=.1;
        row = [zeros(1,i_fix-1),1,zeros(1,inta*intx-i_fix)];
        AT(i_fix,:) = row;

        % Solve linear system
        gg = AT\b;
        g_sum = gg'*ones(inta*intx,1)*da;
        gg = gg./g_sum;
        g = reshape(gg,inta,intx);

        %--------------------------------------------------%
        % Step 1-3 : Checks the stationary equilibrium condition
        %--------------------------------------------------%
        % Update aggregate capital and Check Labor Supply
        S =sum(sum(aa .* g * da)); 
       
        % Excess supply for capital markets
        Ex(iter) = S - K;
       
        if Ex(iter) > critK 
%             disp('Excess Supply')
            rmax = r;
            r = 0.5 * (r + rmin);
        elseif Ex(iter) < -critK
%             disp('Excess Demand')
            rmin = r;
            r = 0.5 * (r + rmax);
        elseif abs(Ex(iter)) < critK
%             disp('Steady state found')
%             disp(r)
            break
        end
        
        % Update Capital and wage based on r
        K = (((alpha) / (r + delta)) ^ (1 / (1 - alpha))) * LAve;
        w = (1 - alpha) * (K ^ alpha) * ((LAve) ^ (-alpha));

        % Initial guess for the value function for the next outer loop : needed???       
        if gamma == 1
            v0(:,1) = log((w * (1 - tau).* x(1) + w * com.* (1 - x(1)) + r.*grida))/rho;
            v0(:,2) = log((w * (1 - tau).* x(2) + w * com.* (1 - x(2)) + r.*grida))/rho;
        else
            v0(:,1) = (w * (1 - tau).* x(1) + w * com.* (1 - x(1)) + r.*grida).^(1-gamma)/(1-gamma)/rho;
            v0(:,2) = (w * (1 - tau).* x(2) + w * com.* (1 - x(2)) + r.*grida).^(1-gamma)/(1-gamma)/rho;
        end
        v = v0;
        
    end % Outer loop
    
    %--------------------------------------------------%
    % Step 2 : Calculates values at the steady state
    %--------------------------------------------------%
    % Values at the deterministic steady state
    rds = r;
    Kds = (((alpha) / (rds + delta)) ^ (1 / (1 - alpha))) * LAve;
    wds = (1 - alpha) * (Kds ^ alpha) * ((LAve) ^ (-alpha));
    cds = c; uds = u; Vds = V; gds = g; 
    pds = wds * (1 - tau).* xx + wds * com.* (1 - xx)  + rds.*aa - c;
    ids = sum((reshape(pds, inta*intx, 1) .* gg)'*da);
end