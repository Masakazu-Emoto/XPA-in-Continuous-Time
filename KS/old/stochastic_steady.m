%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time."
% By Masakazu Emoto and Takeki Sunakawa
% This code calculates stochastic steady state
% This code is refered by Villaverde's code (2019 NBER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ss_mu, Kss, KKss] = stochastic_steady(randZ, muini, Ass, Kdot)
    
    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
    global T N vtime dT
    global a4 x4 K4 Z4
    
    % TFP is zero at stochastic steady state
    if randZ ~= Zmean
        disp('Stop Calculate !')
        stop
    end
    munow = muini;
    
    % Container
    Kss = zeros(N,1); KKss = zeros(N,1);

    for time = 1:N
        % TFP and Aggregate Capital
        Znow = randZ(time); iz = gridsearch(Znow, gridZ);
        if Znow ~= gridZ(iz)
            disp('Stop Calculate !')
            stop
        end
        
        if time == 1
            munext = munow;
        else
            
            % Transition matrix
            Aup = Ass{ss_up,iz}; Adown = Ass{ss_down,iz};
            
            % Calculate wealth distribution
            Mup = (speye(inta*intx) - Aup'*dT);
            Mup = Mup\munow(:);
            Next_up = Mup/(sum(Mup*da));
            Next_up = reshape(Next_up,inta,intx);
            
            Mdown = (speye(inta*intx) - Adown'*dT);
            Mdown = Mdown\munow(:);
            Next_down = Mdown/(sum(Mdown*da));
            Next_down = reshape(Next_down,inta,intx);
            
            % The wealth distribution is calculated by linear interpolation with respect to aggregate capital
            munext = kweight * Next_down + (1 - kweight) * Next_up;
        end
        
        Kss(time) = sum(munext'*grida*da);
        Kss(time) = max([Kss(time) Kmin+0.000001]);
        Kss(time) = min([Kss(time) Kmax-0.000001]);
        
        % Grid search for aggregate capital
        ss_down = floor((Kss(time) - Kmin)/dK) + 1;
        ss_up= ceil((Kss(time) - Kmin)/dK) + 1;
        kweight = (gridK(ss_up) - Kss(time))/dK;
        
        munow = munext;
        
        if time == 1
            Know = sum(muini'*grida*da);
            Knew = Know;
        else
            im = gridsearch(Know, gridK);
            Kweight = (gridK(im + 1) - Know)/dK;
            KK =  (1 - Kweight) * Kdot(im+1, iz) +  Kweight * Kdot(im, iz);

            Knew = Know + KK*dT;
        end
            Know = Knew;
            KKss(time) = Knew;
    end
    
    ss_mu = munext;
end