function [simvalue, sim_mu, simK, simKK] = simulate(randZ, muini, Ass, Kdot)
    
    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
    global T N vtime dT
    
    %%%%%%%%%% 
    % This code is to calculate the stochastic steady state 
    % This code is refered by Villaverde's code (2019 NBER)
    %%%%%%%%%% 
    
    simK = zeros(N,1); simKK = zeros(N,1); 
    
    for time = 1:N
        randZ(time) = max([randZ(time) Zmin+0.000001]);
        randZ(time) = min([randZ(time) Zmax-0.000001]);
        Zdown(time) = floor((randZ(time) - Zmin)/dZ) + 1;
        Zup(time) = ceil((randZ(time) - Zmin)/dZ) + 1;
        zweight(time) = (gridZ(Zup(time)) - randZ(time))/dZ;
    end
        
    munow = muini;

    for time = 1:N
        
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
        
        simK(time) = sum(munext' * grida * da);
        simK(time) = max([simK(time) Kmin+0.000001]);
        simK(time) = min([simK(time) Kmax-0.000001]);
        Kdown(time) = floor((simK(time) - Kmin)/dK) + 1;
        Kup(time) = ceil((simK(time) - Kmin)/dK) + 1;
        kweight = (gridK(Kup(time)) - simK(time))/dK;
        
        munow = munext;
        
        % For calculate Den Haan Error
         if time == 1
            Know = sum(muini'*grida*da);
            Knew = Know;
         else
             
            %             nm = gridsearch(Know, gridK);
            %             nz = gridsearch(randZ(time), gridZ);
            %             Kweight = (gridK(nm + 1) - Know)/dK;
            %             Zweight = (gridZ(nz + 1) - randZ(time))/dZ;
            %             KK =  (1 - Kweight) * (1 - Zweight) * Kdot(nm+1, nz+1) + (1 - Kweight) * Zweight * Kdot(nm+1, nz) + Kweight * (1 - Zweight) * Kdot(nm, nz+1) +  Kweight * Zweight * Kdot(nm, nz);
            
            KKuu = Kdot(KKup(time-1), Zup(time-1));
            KKud = Kdot(KKup(time-1), Zdown(time-1));
            KKdu = Kdot(KKdown(time-1), Zup(time-1));
            KKdd = Kdot(KKdown(time-1), Zdown(time-1));
            
            KK = (1 - Kweight) * (1 - zweight(time - 1)) * KKuu + (1 - Kweight) * zweight(time - 1) * KKud + Kweight * (1 - zweight(time - 1)) * KKdu + Kweight * zweight(time - 1) * KKdd;
            Knew = Know + KK*dT;
         end
            
            simKK(time) = Knew;
            simKK(time) = max([simKK(time) Kmin+0.000001]);
            simKK(time) = min([simKK(time) Kmax-0.000001]);
            
            KKdown(time) = floor((simKK(time) - Kmin)/dK) + 1;
            KKup(time) = ceil((simKK(time) - Kmin)/dK) + 1;
            Kweight = (gridK(KKup(time)) - simKK(time))/dK;
            Know = Knew;
    end
    
    simvalue = zeros(3,1);
    simvalue(1,1) = simK(end);
    simvalue(2,1) = simKK(end);
    simvalue(3,1) = randZ(end);
    
    sim_mu = munext;
end