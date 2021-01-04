function [sim_mu, simK, simKK] = simulate_v1(randZ, muini, Ass, Kdot)
%% simulate_v1.m : This code simulates the complete path of the distribution given the shock, to estimate the forecasting rule
%% INPUTS
% randZ     : the exogenous aggregate shocks
% muini     : the initial distribution from the steady state
% Ass       : the transition matrix obtained in the inner loop
% Kdot      : the forecasting rule
%% OUTPUTS
% sim_mu    : simulated distribution
% simK      : simulated aggregate capital
% simKK     : simulated aggregate capital using the forecasting rule only
% (for the Den Haan error)

    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
    global T N vtime dT
    
    simK = zeros(N,1); simKK = zeros(N,1); 
    
    Zup = zeros(N,1); Zdown = zeros(N,1);
    Kup = zeros(N,1); Kdown = zeros(N,1);
    KKup = zeros(N,1); KKdown = zeros(N,1);
    
    zweight = zeros(N,1); Kweight = zeros(N,1); KKweight = zeros(N,1);
    
    for time = 1:N
        randZ(time)   = max([randZ(time) Zmin+0.000001]);
        randZ(time)   = min([randZ(time) Zmax-0.000001]);
        Zdown(time)   = floor((randZ(time) - Zmin)/dZ) + 1;
        Zup(time)     = ceil((randZ(time) - Zmin)/dZ) + 1;
        zweight(time) = (gridZ(Zup(time)) - randZ(time))/dZ;
    end
        
    %% MAIN LOOP
    munow = muini;

    for time = 1:N
        
        if time == 1
            
            munext = munow;

        else
            
            % Interpolation step to compute next period's distribution
            % The following block computes the transition g(t+1) =
            % (I-A*dt)^(-1)*g(t) for the four closest grid points
            
            % Transition matrix
            Auu = Ass{Kup(time-1), Zup(time-1)}; Aud = Ass{Kup(time-1), Zdown(time-1)};
            Adu = Ass{Kdown(time-1), Zup(time-1)}; Add = Ass{Kdown(time-1), Zdown(time-1)};
            
            Muu = (speye(inta*intx) - Auu'*dT)\munow(:);
            Muu = Muu/(sum(Muu*da));
            Nextuu = reshape(Muu,inta,intx);
            
            Mud = (speye(inta*intx) - Aud'*dT)\munow(:);
            Mud = Mud/(sum(Mud*da));
            Nextud = reshape(Mud,inta,intx);
            
            Mdu = (speye(inta*intx) - Adu'*dT)\munow(:);
            Mdu = Mdu/(sum(Mdu*da));
            Nextdu = reshape(Mdu,inta,intx);
            
            Mdd = (speye(inta*intx) - Add'*dT)\munow(:);
            Mdd = Mdd/(sum(Mdd*da));
            Nextdd = reshape(Mdd,inta,intx);
            
            % Then we compute the wealth distribution in the next period, averaging these four results by linear interpolation with respect to aggregate capital and aggregate productivity
            munext = (1 - Kweight) * (1 - zweight(time - 1)) * Nextuu + (1 - Kweight) * zweight(time - 1) * Nextud + Kweight * (1 - zweight(time - 1)) * Nextdu + Kweight * zweight(time - 1) * Nextdd;
            
        end
        
        simK(time) = sum(munext' * grida * da);

        simK(time) = max([simK(time) Kmin+0.000001]);
        simK(time) = min([simK(time) Kmax-0.000001]);
        Kdown(time) = floor((simK(time) - Kmin)/dK) + 1;
        Kup(time) = ceil((simK(time) - Kmin)/dK) + 1;
        Kweight = (gridK(Kup(time)) - simK(time))/dK;
        
        munow = munext; % update the distribution
        
        %% For calculating the Den Haan Error
        % simKK is simulated using the forecasting rule only
        if time == 1

            Know = sum(muini'*grida*da);
            Knew = Know;

        else
            
            % Interpolation using the four closest grid points
            KKuu = Kdot(KKup(time-1), Zup(time-1));
            KKud = Kdot(KKup(time-1), Zdown(time-1));
            KKdu = Kdot(KKdown(time-1), Zup(time-1));
            KKdd = Kdot(KKdown(time-1), Zdown(time-1));
            
            % KKweight is based on simKK(time-1)
            KK = (1 - KKweight) * (1 - zweight(time - 1)) * KKuu + (1 - KKweight) * zweight(time - 1) * KKud + KKweight * (1 - zweight(time - 1)) * KKdu + KKweight * zweight(time - 1) * KKdd;
            Knew = Know + KK*dT;
        
        end
            
        simKK(time) = Knew;

        simKK(time) = max([simKK(time) Kmin+0.000001]);
        simKK(time) = min([simKK(time) Kmax-0.000001]);
        KKdown(time) = floor((simKK(time) - Kmin)/dK) + 1;
        KKup(time) = ceil((simKK(time) - Kmin)/dK) + 1;
        KKweight = (gridK(KKup(time)) - simKK(time))/dK;
        
        Know = Knew;
    
    end
    
    sim_mu = munext;
    
end