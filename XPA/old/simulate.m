%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time."
% By Masakazu Emoto and Takeki Sunakawa
% This code simulates and derives the path of aggregate capital
% This code is refered by Villaverde's code (2019 NBER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [sim_mu, simK, simKK] = simulate(randZ, muini, Ass, Kdot)

    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
    global T N vtime dT
    
    % Container
    simK = zeros(N,1); simKK = zeros(N,1); 
    
    Zup = zeros(N,1); Zdown = zeros(N,1);
    Kup = zeros(N,1); Kdown = zeros(N,1);
    KKup = zeros(N,1); KKdown = zeros(N,1);
    
    zweight = zeros(N,1); Kweight = zeros(N,1); KKweight = zeros(N,1);
    
    % Grid search for aggregate uncertainty
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
            
            % Transition matrix
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
            
            % The wealth distribution is calculated by linear interpolation with respect to aggregate capital and aggregate productivity.
            munext = (1 - Kweight) * (1 - zweight(time - 1)) * Nextuu + (1 - Kweight) * zweight(time - 1) * Nextud + Kweight * (1 - zweight(time - 1)) * Nextdu + Kweight * zweight(time - 1) * Nextdd;
        end
        
        % simK is simulated results using the dynamics of aggregate capital and the HJB equation
        simK(time) = sum(munext' * grida * da);
        simK(time) = max([simK(time) Kmin+0.000001]);
        simK(time) = min([simK(time) Kmax-0.000001]);
        
        % Grid search for aggregate capital
        Kdown(time) = floor((simK(time) - Kmin)/dK) + 1;
        Kup(time) = ceil((simK(time) - Kmin)/dK) + 1;
        Kweight = (gridK(Kup(time)) - simK(time))/dK;
        
        munow = munext;
        
        % For calculate Den Haan Error
         if time == 1
            Know = sum(muini'*grida*da);
            Knew = Know;
         else
            % Calculate law of motion
            KKuu = Kdot(KKup(time-1), Zup(time-1));
            KKud = Kdot(KKup(time-1), Zdown(time-1));
            KKdu = Kdot(KKdown(time-1), Zup(time-1));
            KKdd = Kdot(KKdown(time-1), Zdown(time-1));
            
            % The law of motion is calculated by linear interpolation with respect to aggregate capital and aggregate productivity.
            KK = (1 - KKweight) * (1 - zweight(time - 1)) * KKuu + (1 - KKweight) * zweight(time - 1) * KKud + KKweight * (1 - zweight(time - 1)) * KKdu + KKweight * zweight(time - 1) * KKdd;
            Knew = Know + KK*dT;
         end
            
        % simKK is simulated results using the dynamics of aggregate capital
        simKK(time) = Knew;
        simKK(time) = max([simKK(time) Kmin+0.000001]);
        simKK(time) = min([simKK(time) Kmax-0.000001]);
        
        % Grid search for aggregate capital
        KKdown(time) = floor((simKK(time) - Kmin)/dK) + 1;
        KKup(time) = ceil((simKK(time) - Kmin)/dK) + 1;
        KKweight = (gridK(KKup(time)) - simKK(time))/dK;
        Know = Knew;
    end
    
    sim_mu = munext;
end