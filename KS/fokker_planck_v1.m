function [Ksim, Zsim, Kdown, Kup, Zdown, Zup] = fokker_planck_v1(Zshocks, muini, Ass)
%% fokker_planck_v1.m : This code simulates the complete path of the distribution given the shock, to estimate the forecasting rule
%% INPUTS
% Zshocks   : the exogenous aggregate shocks
% muini     : the initial distribution from the steady state
% Ass       : the transition matrix obtained in the inner loop
%% OUTPUTS
% Ksim      : simulated aggregate capital
% Zsim      : simulated aggregate productivity
% Kdown     : the downward closest point of Kgrid
% Kup       : the upward closest point of Kgrid
% Zdown     : the downward closest point of Zgrid
% Zup       : the upward closest point of Zgrid

%% NOTE: This code is based on b5_KFE.m written by FVHN (2018).

    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK intKK gridK dK dKK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
    global T N Stime Dtime vtime dT
 
    %% This part is to be moved to the main file
    Zsim    = zeros(Stime,1);
    Zup     = zeros(Stime,1);
    Zdown   = zeros(Stime,1);
    zweight = zeros(Stime,1);

    % The path of aggregate productvity
    for time = 1:Stime-1
        if time == 1
            Zsim(time+1) = mu *dT * Zmean + (1 - mu * dT) * Zmean + sigma * Zshocks(time) * sqrt(dT);
        else
            Zsim(time+1) = mu *dT * Zmean + (1 - mu * dT) * Zsim(time) + sigma * Zshocks(time) * sqrt(dT);
        end
    end
    
    for time = 1:Stime
        Zsim(time)    =   max([Zsim(time) Zmin+0.000001]);
        Zsim(time)    =   min([Zsim(time) Zmax-0.000001]);
        Zdown(time)   = floor((Zsim(time)-Zmin)/dZ)+1;
        Zup(time)     =  ceil((Zsim(time)-Zmin)/dZ)+1;
        zweight(time) = (gridZ(Zup(time))-Zsim(time))/dZ;    % weight of ZposD
    end

    %% MAIN LOOP
    munow = muini;

    Ksim   = zeros(Stime,1);
    Kup    = zeros(Stime,1); % KposU
    Kdown  = zeros(Stime,1); % KposD
%     KKup   = zeros(Stime,1); % KKposU
%     KKdown = zeros(Stime,1); % KKposD

    for time = 1:Stime
        
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
            munext = (1 - kweight) * (1 - zweight(time - 1)) * Nextuu + (1 - kweight) * zweight(time - 1) * Nextud + kweight * (1 - zweight(time - 1)) * Nextdu + kweight * zweight(time - 1) * Nextdd;
            
        end
        
        Ksim(time) = sum(munext' * grida * da);
        
        Ksim(time) = max([Ksim(time) Kmin+0.000001]);
        Ksim(time) = min([Ksim(time) Kmax-0.000001]);
        Kdown(time) = floor((Ksim(time) - Kmin)/dK) + 1;
        Kup(time) = ceil((Ksim(time) - Kmin)/dK) + 1;
        kweight = (gridK(Kup(time)) - Ksim(time))/dK;
        
        munow = munext; % update the distribution
        
    end
    
end