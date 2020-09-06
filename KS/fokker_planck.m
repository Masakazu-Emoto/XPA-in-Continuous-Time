function [Ksim, Zsim, Zdown, Zup, Kdown, Kup] = fokker_planck(Zshocks, muini, Ass)
    
    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK intKK gridK dK dKK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
    global T N Stime Dtime vtime dT
    
    %%%%%%%%%% 
    % This code is to calculate the stochastic steady state 
    % This code is refered by Villaverde's code (2019 NBER)
    %%%%%%%%%% 
    
    % Aggregate Productivity
    Zsim    = zeros(Stime,1);
    Zup   = zeros(Stime,1);
    Zdown   = zeros(Stime,1);
    zweight      = zeros(Stime,1);

    for time = 1:Stime-1
        if time == 1
            Zsim(time+1) = mu *dT * Zmean + (1 - mu * dT) * Zmean + sigma * Zshocks(time) * sqrt(dT);
        else
            Zsim(time+1) = mu *dT * Zmean + (1 - mu * dT) * Zsim(time) + sigma * Zshocks(time) * sqrt(dT);
        end
    end
    
    for time = 1:Stime
        Zsim(time)  =   max([Zsim(time) Zmin+0.000001]);
        Zsim(time)  =   min([Zsim(time) Zmax-0.000001]);
        Zdown(time) = floor((Zsim(time)-Zmin)/dZ)+1;
        Zup(time) =  ceil((Zsim(time)-Zmin)/dZ)+1;
        zweight(time)    = (gridZ(Zup(time))-Zsim(time))/dZ;    % weight of ZposD
    end

    munow = muini;
    % Aggregate Productivity
    Ksim = zeros(Stime,1);
    Kup   = zeros(Stime,1);
    Kdown   = zeros(Stime,1);
    KKup = zeros(Stime,1);
    KKdown = zeros(Stime,1);

    for time = 1:Stime
        
        if time == 1
            munext = munow;
        else
            
            % Calculate individual consumption, saving and coefficient wealth distribution)};
            Auu = Ass{Kup(time-1), Zup(time-1)}; Aud = Ass{Kup(time-1), Zdown(time-1)};
            Adu = Ass{Kdown(time-1), Zup(time-1)}; Add = Ass{Kdown(time-1), Zdown(time-1)};
            
            % Calculate wealth distribution
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
            
            munext = (1 - kweight) * (1 - zweight(time - 1)) * Nextuu + (1 - kweight) * zweight(time - 1) * Nextud + kweight * (1 - zweight(time - 1)) * Nextdu + kweight * zweight(time - 1) * Nextdd;
            
        end
        
%         Ksim(time) = sum(munext' * grida * da);
        Ksim(time) = sum(sum(munext.*aa*da));
        Ksim(time) = max([Ksim(time) Kmin+0.000001]);
        Ksim(time) = min([Ksim(time) Kmax-0.000001]);
        Kdown(time) = floor((Ksim(time) - Kmin)/dK) + 1;
        Kup(time) = ceil((Ksim(time) - Kmin)/dK) + 1;
        kweight = (gridK(Kup(time)) - Ksim(time))/dK;
        
        munow = munext;
    end
end