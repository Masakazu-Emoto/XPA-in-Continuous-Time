function [A1, A1tilde, A3, vss, cs, ps, zx, zy, zz] = inner_org(Kdot, vss, r, w)
%% inner_org.m : This code solves the Hamilton-Jacobi-Bellman equation with aggregate uncertainty by the finite differential method
%% INPUTS
% Kdot      : the forecasting rule
% vss       : the initial value function from the steady state
% r         : the real interest rate
% w         : the wage rate
%% OUTPUTS
% A1    : A_lm
% A1tilde : A_lm, excluding the effect of aggregate variables, to be used when we solve the KF equation
% A3    : the transition matrix, A3
% cs    : The policy function for consumption
% ps    : The policy function for savings
% zx    : Forward difference of V in terms of Z (to be removed)
% zy    : Central difference of V in terms of Z
% zz    : Backward difference of V in terms of Z

%% NOTE: This code is based on b3_HJB.m written by FVHN (2018).
% However, we extend their original code in the following two dimensions:
% (1) We use the upwind scheme not only individual wealth, a, but also K and Z.
% (2) We exclude the direct effect of aggregate variables K and Z on the matrix
% A_lm in their note when we solve the KF equation (there is the indirect effect through r and w).

    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
    global quada quadx quadK quadZ
    
    % Finite difference approximation of the partial derivatives
    % 4-dimentional matrices, the ones we want to reset every iteration
    Vsaf = zeros(inta,intx,intK,intZ); Vsab = zeros(inta,intx,intK,intZ);
    dVdK = zeros(inta,intx,intK,intZ); dVdZ = zeros(inta,intx,intK,intZ); 
    dVddZ = zeros(inta,intx,intK,intZ);

    % for consistency (to be moved to main file)
    zx = -min(zmu,0)/dZ + zsigma/(2*ddZ);
    zy = min(zmu,0)/dZ - max(zmu, 0)/dZ - zsigma/ddZ;
    zz = max(zmu,0)/dZ + zsigma/(2*ddZ);

    % Collections of empty sparse matrices (cell arrays)    
    % A1 could also be called A_lm, A2 could also be called A_m, A3 could also be called A    
    % Ass = A1 and WW = A3 in FVHN (to be fixed)
    for ik=1:intK
        for iz=1:intZ
            A1{ik,iz}=sparse(inta*intx, inta*intx);
            A1tilde{ik,iz}=sparse(inta*intx, inta*intx);
        end
    end

    for iz=1:intZ
        A2{iz}=sparse(inta*intx*intK, inta*intx*intK);
    end

    A3 = sparse(inta*intx*intK*intZ, inta*intx*intK*intZ);

    % Law of motion matirix
    PLM = zeros(inta,intx,intK,intZ);
    for ik = 1:intK
        for iz = 1:intZ
            PLM(:,:,ik,iz) = Kdot(ik,iz);
        end
    end
    
    % TS: BAD NOTATION (WHY "ss"??????)
    Vss = vss;

    % MAIN LOOP
    for n = 1:maxit

        Vss = 0.5 * vss + (1 - 0.5) * Vss;
        % Forward Difference : Individual wealth
        Vsaf(1:inta-1,:,:,:) = (Vss(2:inta,:,:,:)-Vss(1:inta-1,:,:,:))/da;
        Vsaf(inta,:,:,:) = (w(inta,:,:,:) * (1 - tau).*quadx(inta,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(inta,:,:,:))+ r(inta,:,:,:).*amax).^(-gamma);

        % Backward Difference : Individual wealth
        Vsab(2:inta,:,:,:) = (Vss(2:inta,:,:,:)-Vss(1:inta-1,:,:,:))/da; %w(inta,:,:,:) * (1 - tau).*quadx(inta,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(inta,:,:,:))
        Vsab(1,:,:,:) = (w(1,:,:,:) * (1 - tau).*quadx(1,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(1,:,:,:)) + r(1,:,:,:).*amin).^(-gamma);

        % Finite Difference : Aggregate K and Z
        dVdK(:,:,1:intK-1,:)  = (Vss(:,:,2:intK,:)-Vss(:,:,1:intK-1,:))/dK;
        dVdZ(:,:,:,1:intZ-1)  = (Vss(:,:,:,2:intZ)-Vss(:,:,:,1:intZ-1))/dZ;
        dVddZ(:,:,:,2:intZ-1) = (Vss(:,:,:,3:intZ)-2*Vss(:,:,:,2:intZ-1)+Vss(:,:,:,1:intZ-2))/(dZ^2);
        
        % Consumption and savings with forward difference
        cf = (Vsaf).^(-1/gamma); 
        sf = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cf;

        % Consumption and savings with backward difference
        cb = (Vsab).^(-1/gamma);
        sb = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cb;

        % Consumption and derivative of value function at steady state
        c0 = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada;
        Va0 = c0.^(-gamma);

        % dV_upwind makes a choice of forward or backward differences based on
        % The sign of the drift for individual wealth
        Iaf = (sf > 0);         % Positive drift --> Forward difference
        Iab = (sb < 0);         % Negative drift --> Backward difference
        Ia0 = (1 - Iaf - Iab);  % Drift is zero --> At steady state

        Va_Upwind = Vsaf.*Iaf + Vsab.*Iab + Va0.*Ia0; % important to include third term
        cs = (Va_Upwind).^(-1/gamma); ps = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cs; 
        if gamma == 1
            uss = log(cs);
        else
            uss = (cs.^(1-gamma))/(1-gamma);
        end
        
        % -------------------------------------------------- %    
        % FVHN Algorithm
        % -------------------------------------------------- % 
        % Construct sparse matrix
        elem_a       = (-sb.*Iab)/da;
        elem_e       = ( sf.*Iaf)/da;
        elem_b       = -elem_a-elem_e - PLM/dK - mu*(Zmean-quadZ)/dZ - ((sigma/dZ)^2);
        elem_btilde  = -elem_a-elem_e; % for A1tilde matrix
        elem_r       = ((sigma/dZ)^2)/2;                % this one is a scalar
        elem_x       = mu*(Zmean-gridZ)/dZ + elem_r;    % this one is a vector
        
        % Fill A1 collection
        for ik=1:intK
            for iz=1:intZ
                A11          = spdiags(elem_b(:,1,ik,iz),0,inta,inta) + spdiags(elem_a(2:inta,1,ik,iz),-1,inta,inta) + spdiags([0;elem_e(1:inta-1,1,ik,iz)],1,inta,inta);
                A12          = spdiags(elem_b(:,2,ik,iz),0,inta,inta) + spdiags(elem_a(2:inta,2,ik,iz),-1,inta,inta) + spdiags([0;elem_e(1:inta-1,2,ik,iz)],1,inta,inta);
                A1{ik,iz}    = [A11,sparse(inta, inta); sparse(inta, inta),A12] + Aswitch;
                A11tilde     = spdiags(elem_btilde(:,1,ik,iz),0,inta,inta) + spdiags(elem_a(2:inta,1,ik,iz),-1,inta,inta) + spdiags([0;elem_e(1:inta-1,1,ik,iz)],1,inta,inta);
                A12tilde     = spdiags(elem_btilde(:,2,ik,iz),0,inta,inta) + spdiags(elem_a(2:inta,2,ik,iz),-1,inta,inta) + spdiags([0;elem_e(1:inta-1,2,ik,iz)],1,inta,inta);
                A1tilde{ik,iz} = [A11tilde,sparse(inta, inta); sparse(inta, inta),A12tilde] + Aswitch;
            end
        end
        
        % Fill A2 collection
        jumper =inta*intx;
        for iz = 1:intZ
            for ik = 1:intK
                A2{iz}((ik - 1) * jumper + 1:ik * jumper,(ik - 1) * jumper + 1:ik * jumper) = A1{ik,iz};
            end
        end
        for iz = 1:intZ
            for ik = 1:intK - 1
                A2{iz}((ik - 1)*jumper+1:ik*jumper, ik*jumper+1:(ik + 1)*jumper) = (Kdot(ik,iz)/dK) * speye(inta*intx,inta*intx);
            end
            A2{iz}((intK - 1)*jumper+1:intK*jumper, (intK - 1)*jumper + 1:intK*jumper)=A2{iz}((intK-1)*jumper+1:intK*jumper,(intK-1)*jumper+1:intK*jumper)+(Kdot(intK,iz)/dK)*speye(inta*intx,inta*intx);
        end
        
        % Fill A3
        jumper = inta*intx*intK;
        for iz=1:intZ
            A3((iz-1)*jumper+1:iz*jumper,(iz-1)*jumper+1:iz*jumper) = A2{iz};
        end
        for iz=1:intZ-1
            A3(iz*jumper+1:(iz+1)*jumper,(iz-1)*jumper+1:iz*jumper) = elem_r*speye(inta*intx*intK);
        end
        for iz=1:intZ-1
            A3((iz-1)*jumper+1:iz*jumper,iz*jumper+1:(iz+1)*jumper) = elem_x(iz)*speye(inta*intx*intK);
        end
        A3((intZ-1)*jumper+1:intZ*jumper,(intZ-1)*jumper+1:intZ*jumper) = A3((intZ-1)*jumper+1:intZ*jumper,(intZ-1)*jumper+1:intZ*jumper) + elem_x(intZ)*speye(inta*intx*intK, inta*intx*intK);
        A3(1:jumper,1:jumper) = A3(1:jumper,1:jumper) + elem_r * speye(inta*intx*intK);
        
        % Check if the transition matrix is real
        if ~isreal(A3)
            disp('')
            disp('The transition matrix is not real')
            disp('')
            break
        end
        
        % Calculate BB and bb to find Vnew
        BB = (rho + 1/Delta) * speye(inta*intx*intK*intZ) - A3;
        u_stacked = uss(:); V_stacked = Vss(:);

        bb = u_stacked + V_stacked/Delta;

        Vnew_stacked = BB\bb; % SOLVE SYSTEM OF EQUATIONS
        
        vss = reshape(Vnew_stacked,size(Vss));
        Vchange = Vnew_stacked - V_stacked;
        
        dist = max(abs(Vchange));
        
        % Check convergence
        if dist < crit
%             disp('Value Function Converged')
%             disp('Difference = ')
%             disp(dist); 
%             disp('Iteration = ')
%             disp(n); 
            
            % If the value function is converged, HJB_check is as close to zero as possible.
            % FVHN case
            HJB_check = zeros(inta, intx, intK, intZ);
            if gamma == 1
                HJB_check(:,1,:,:) = -rho*vss(:,1,:,:) + log(cs(:,1,:,:)) + (w(:,1,:,:)* (1 - tau).*quadx(:,1,:,:) + w(:,1,:,:)* com.*(1 - quadx(:,1,:,:)) + r(:,1,:,:).*quada(:,1,:,:)-cs(:,1,:,:)).*Va_Upwind(:,1,:,:) + la(1)*(vss(:,2,:,:)-vss(:,1,:,:)) + PLM(:,1,:,:).*dVdK(:,1,:,:) + mu*(Zmean-quadZ(:,1,:,:)).*dVdZ(:,1,:,:) + ((sigma^2)/2)*dVddZ(:,1,:,:);
                HJB_check(:,2,:,:) = -rho*vss(:,2,:,:) + log(cs(:,2,:,:)) + (w(:,2,:,:)* (1 - tau).*quadx(:,2,:,:) + w(:,2,:,:)* com.*(1 - quadx(:,2,:,:)) + r(:,2,:,:).*quada(:,2,:,:)-cs(:,2,:,:)).*Va_Upwind(:,2,:,:) + la(2)*(vss(:,1,:,:)-vss(:,2,:,:)) + PLM(:,2,:,:).*dVdK(:,2,:,:) + mu*(Zmean-quadZ(:,2,:,:)).*dVdZ(:,2,:,:) + ((sigma^2)/2)*dVddZ(:,2,:,:);
            else
                HJB_check(:,1,:,:) = -rho*vss(:,1,:,:) + (cs(:,1,:,:).^(1-gamma))/(1-gamma) + (w(:,1,:,:)* (1 - tau).*quadx(:,1,:,:) + w(:,1,:,:)* com.*(1 - quadx(:,1,:,:)) + r(:,1,:,:).*quada(:,1,:,:)-cs(:,1,:,:)).*Va_Upwind(:,1,:,:) + la(1)*(vss(:,2,:,:)-vss(:,1,:,:)) + PLM(:,1,:,:).*dVdK(:,1,:,:) + mu*(Zmean-quadZ(:,1,:,:)).*dVdZ(:,1,:,:) + ((sigma^2)/2)*dVddZ(:,1,:,:);
                HJB_check(:,2,:,:) = -rho*vss(:,2,:,:) + (cs(:,2,:,:).^(1-gamma))/(1-gamma) + (w(:,2,:,:)* (1 - tau).*quadx(:,2,:,:) + w(:,2,:,:)* com.*(1 - quadx(:,2,:,:)) + r(:,2,:,:).*quada(:,2,:,:)-cs(:,2,:,:)).*Va_Upwind(:,2,:,:) + la(2)*(vss(:,1,:,:)-vss(:,2,:,:)) + PLM(:,2,:,:).*dVdK(:,2,:,:) + mu*(Zmean-quadZ(:,2,:,:)).*dVdZ(:,2,:,:) + ((sigma^2)/2)*dVddZ(:,2,:,:);
            end
                       
            HJB_check_stacked = zeros(inta*intx*intK*intZ,1);
            for iZ=2:intZ-1 % reshape could be faster, but this is clearer
                for iK=2:intK-1
                    for ix=1:intx
                        HJB_check_stacked(2+(ix-1)*inta+(iK-1)*inta*intx+(iZ-1)*inta*intx*intK:(inta-1)+(ix-1)*inta+(iK-1)*inta*intx+(iZ-1)*inta*intx*intK,1)=HJB_check(2:inta-1,ix,iK,iZ);
                    end
                end
            end

%             disp('HJB_check:')
%             disp(max(max(abs(HJB_check_stacked))))
            break
        end
    end

    % Maybe redundant - cs and ps are already obtained
%     % -------------------------------------------------- %
%     % Step 6 : Calculate Consumption function and Policy function
%     % -------------------------------------------------- %
%     % Forward difference : Indivudual wealth
%     Vsaf(1:inta-1,:,:,:) = (Vss(2:inta,:,:,:)-Vss(1:inta-1,:,:,:))/da;
%     Vsaf(inta,:,:,:) = (w(inta,:,:,:) * (1 - tau).*quadx(inta,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(inta,:,:,:))+ r(inta,:,:,:).*amax).^(-gamma);
% 
%     % Backward difference : Individual wealth
%     Vsab(2:inta,:,:,:) = (Vss(2:inta,:,:,:)-Vss(1:inta-1,:,:,:))/da;
%     Vsab(1,:,:,:) = (w(1,:,:,:) * (1 - tau).*quadx(1,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(1,:,:,:)) + r(1,:,:,:).*amin).^(-gamma);
%     
%     % Consumption and savings with forward difference
%     cf = (Vsaf).^(-1/gamma); 
%     sf = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cf;
% 
%     % Consumption and savings with backward difference
%     cb = (Vsab).^(-1/gamma);
%     sb = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cb;
% 
%     % Consumption and derivative of value function at steady state
%     c0 = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada;
%     Va0 = c0.^(-gamma);
% 
%     % dV_upwind makes a choice of forward or backward differences based on
%     % The sign of the drift for individual wealth
%     Iaf = (sf > 0);             % Positive drift --> Forward difference
%     Iab = (sb < 0);           % Negative drift --> Backward difference
%     Ia0 = (1 - Iaf - Iab);   % % Drift is zero --> At steady state
% 
%     Va_Upwind = Vsaf.*Iaf + Vsab.*Iab + Va0.*Ia0;
%     cs = (Va_Upwind).^(-1/gamma); ps = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cs;
    
end