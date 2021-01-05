function [Ksim, KKsim] = simulate_v2(Zsim, Zdown, Zup, zweight, muini, Ass, Kdot)
%% simulate_v2.m : This code simulates the complete path of the distribution given the shock, to estimate the forecasting rule
%
%% INPUTS
%    Zsim       : simulated aggregate productivity
%    Zdown      : the downward closest point of Zgrid
%    Zup        : the upward closest point of Zgrid
%    zweight    : weight based on a value of Z and Zgrid
%    muini      : the initial distribution from the steady state
%    Ass        : the transition matrix obtained in the inner loop
%    Kdot       : the forecasting rule
%
%% OUTPUTS
%    Ksim       : simulated aggregate capital
%    KKsim      : simulated aggregate capital using the forecasting rule only (for the Den Haan error)
%
%% NOTE: This code is based on b5_KFE.m written by FVHN.

    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
%    global T N vtime dT
    global dT
    
    N = size(Zsim,1);
        
    %% MAIN LOOP
    munow = muini;

    Ksim  = zeros(N,1); 
    KKsim = zeros(N,1); 
    Kup   = zeros(N,1); 
    Kdown = zeros(N,1);
    KKup   = zeros(N,1); 
    KKdown = zeros(N,1);
    
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
        
        Ksim(time)  = sum(munext' * grida * da);

        Ksim(time)  = max([Ksim(time) Kmin+0.000001]);
        Ksim(time)  = min([Ksim(time) Kmax-0.000001]);
        Kdown(time) = floor((Ksim(time) - Kmin)/dK) + 1;
        Kup(time)   = ceil((Ksim(time) - Kmin)/dK) + 1;
        Kweight     = (gridK(Kup(time)) - Ksim(time))/dK;
        
        munow = munext; % update the distribution
        
        %% For calculating the Den Haan Error
        if (nargin>6)
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

            KKsim(time)  = Knew;

            KKsim(time)  = max([KKsim(time) Kmin+0.000001]);
            KKsim(time)  = min([KKsim(time) Kmax-0.000001]);
            KKdown(time) = floor((KKsim(time) - Kmin)/dK) + 1;
            KKup(time)   = ceil((KKsim(time) - Kmin)/dK) + 1;
            KKweight     = (gridK(KKup(time)) - KKsim(time))/dK;

            Know = Knew;
        
        end
    
    end
        
end