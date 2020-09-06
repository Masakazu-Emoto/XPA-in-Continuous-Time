% function [vss, pss, Zup, Zcenter, Zdown, Kup, Kcenter, Kdown] = inner(Kdot)
function [Ass, WW, vss, cs, ps] = inner(Kdot, vss, iteration, r, w)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Algorithm : Inner loop by XPA in continuous time
    % Phase 1 : Construct sparse matrix for aggregate shock
    % Phase 2 : Construct sparse matrix for aggregate wealth
    % Phase 3 : Initial guess for some variables
    % Phase 4 : Calculating value function and policy function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global gamma rho alpha delta la intx x mu sigma com tau LAve 
    global maxit maxitK crit critK Delta damp
    global inta amin amax grida da aa aaa xx xxx Aswitch
    global Kmax Kmin intK gridK dK
    global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ
    global  quada quadx quadK quadZ
    
    %Finite difference approximation of the partial derivatives
    Vsaf = zeros(inta,intx,intK,intZ); Vsab = zeros(inta,intx,intK,intZ);
    dVdK   = zeros(inta,intx,intK,intZ); dVdZ   = zeros(inta,intx,intK,intZ); dVddZ = zeros(inta,intx,intK,intZ);
    
    % -------------------------------------------------- %
    % Phase 1 : Construct sparse matrix for aggregate shock
    zx = -min(zmu,0)/dZ + zsigma/(2*ddZ);
    zy = min(zmu,0)/dZ - max(zmu, 0)/dZ - zsigma/ddZ;
    zz = max(zmu,0)/dZ + zsigma/(2*ddZ);
    
    % This will be the upperdiagonal of the matrix
    Zup = cell(intK,intZ); Zcenter = cell(intK,intZ); Zdown = cell(intK,intZ);
    for iz = 1:intZ
        for ik = 1:intK
            if iz == 1
                Zup{ik,iz} = zz(iz,1)*speye(inta*intx, inta*intx); 
                Zcenter{ik,iz} = (zy(iz,1) + zx(iz,1))*speye(inta*intx, inta*intx);
                Zdown{ik,iz} = sparse(inta*intx, inta*intx);
            elseif iz == intZ
                Zup{ik,iz} = sparse(inta*intx, inta*intx);
                Zcenter{ik,iz} = (zy(iz,1) + zz(iz,1))*speye(inta*intx, inta*intx);
                Zdown{ik,iz} = zx(iz,1)*speye(inta*intx, inta*intx);
            else    
                Zup{ik,iz} = zz(iz,1)*speye(inta*intx, inta*intx);
                Zcenter{ik,iz} = zy(iz,1)*speye(inta*intx, inta*intx);
                Zdown{ik,iz} = zx(iz,1)*speye(inta*intx, inta*intx);
            end
        end
    end
    
    %For calculate
    ZZup = reshape(Zup,intK*intZ,1);
    ZZdown = reshape(Zdown,intK*intZ,1);
    
    % -------------------------------------------------- %
    % Phase 2 : Construct sparse matrix for aggregate wealth
    kx = -min(Kdot,0)/dK;
    ky = min(Kdot,0)/dK - max(Kdot,0)/dK;
    kz = max(Kdot,0)/dK;
    
    % This will be the upperdiagonal of the matrix 
    Kup = cell(intK,intZ); Kcenter = cell(intK,intZ); Kdown = cell(intK,intZ);
    for iz = 1:intZ
        for ik = 1:intK
            if ik == 1
                Kup{ik,iz} = kz(ik,iz)*speye(inta*intx,inta*intx);
                Kcenter{ik,iz} = (ky(ik,iz) + kx(ik,iz))*speye(inta*intx,inta*intx);
                Kdown{ik,iz} = sparse(inta*intx,inta*intx);
            elseif ik == intK
                Kup{ik,iz} = sparse(inta*intx,inta*intx);
                Kcenter{ik,iz} = (ky(ik,iz) + kz(ik,iz))*speye(inta*intx,inta*intx);
                Kdown{ik,iz} = kx(ik,iz)*speye(inta*intx,inta*intx);
            else
                Kup{ik,iz} = kz(ik,iz)*speye(inta*intx,inta*intx);
                Kcenter{ik,iz} = ky(ik,iz)*speye(inta*intx,inta*intx);
                Kdown{ik,iz} = kx(ik,iz)*speye(inta*intx,inta*intx);
            end
        end
    end
    
    KKup = reshape(Kup,intK*intZ,1);
    KKdown = reshape(Kdown,intK*intZ,1);
    
    % Villaverde Algorithm
%         for ik=1:intK
%             for iz=1:intZ
%                 Ass{ik,iz}=sparse(inta*intx, inta*intx);
%             end
%         end
% 
%         for iz=1:intZ
%             A2{iz}=sparse(inta*intx*intK, inta*intx*intK);
%         end
% 
%         WW = sparse(inta*intx*intK*intZ, inta*intx*intK*intZ);
%         
        PLM = zeros(inta,intx,intK,intZ);
        for ik = 1:intK
            for iz = 1:intZ
                PLM(:,:,ik,iz) = Kdot(ik,iz);
            end
        end
    
    % -------------------------------------------------- %
    % Phase 3 : Initial guess for some variable
    % Initial value function
    Vss = vss;
    
    % -------------------------------------------------- %    
    % Phase 4 : Calculating value function and policy function
    for n = 1:maxit

        Vss = 0.5 * vss + (1 - 0.5) * Vss;
        % Forward Difference : Indivudual wealth
        Vsaf(1:inta-1,:,:,:) = (Vss(2:inta,:,:,:)-Vss(1:inta-1,:,:,:))/da;
        Vsaf(inta,:,:,:) = (w(inta,:,:,:) * (1 - tau).*quadx(inta,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(inta,:,:,:))+ r(inta,:,:,:).*amax).^(-gamma);
%         Vsaf(inta,:,:,:) = (w(inta,:,:,:).*quadx(inta,:,:,:) + r(inta,:,:,:).*amax).^(-gamma);

        % Backward Difference : Indivisual wealth
        Vsab(2:inta,:,:,:) = (Vss(2:inta,:,:,:)-Vss(1:inta-1,:,:,:))/da; %w(inta,:,:,:) * (1 - tau).*quadx(inta,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(inta,:,:,:))
        Vsab(1,:,:,:) = (w(1,:,:,:) * (1 - tau).*quadx(1,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(1,:,:,:)) + r(1,:,:,:).*amin).^(-gamma);
%         Vsab(1,:,:,:) = (w(1,:,:,:).*quadx(1,:,:,:) + r(1,:,:,:).*amin).^(-gamma);

        % Consumption and savings with forward difference
        cf = (Vsaf).^(-1/gamma); 
        sf = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cf;
%         sf = w.*quadx + r.*quada - cf;

        % Consumption and savings with backward difference
        cb = (Vsab).^(-1/gamma);
        sb = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cb;
%         sb = w.*quadx + r.*quada - cb;

        % Consumption and derivative of value function at steady state
        c0 = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada;
        Va0 = c0.^(-gamma);
%         c0 = w.*quadx + r.*quada;

        % dV_upwind makes a choice of forward or backward differences based on
        % The sign of the drift for individual wealth
        Iaf = (sf > 0);             % Positive drift --> Forward difference
        Iab = (sb < 0);           % Negative drift --> Backward difference
        Ia0 = (1 - Iaf - Iab); % At steady state

        Va_Upwind = Vsaf.*Iaf + Vsab.*Iab + Va0.*Ia0; %important to include third term
        cs = (Va_Upwind).^(-1/gamma); ps = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cs; 
        if gamma == 1
            uss = log(cs);
        else
            uss = (cs.^(1-gamma))/(1-gamma);
        end
%     ps = w.*quadx + r.*quada - cs;

        % Villaverde Algorithm
%         elem_a       = (-sb.*Iab)/da;
%         elem_e       = ( sf.*Iaf)/da;
%         elem_b       = -elem_a-elem_e - PLM/dK - mu*(Zmean-quadZ)/dZ - ((sigma/dZ)^2);
%         elem_r       = ((sigma/dZ)^2)/2;                       % this one is a scalar
%         elem_x       = mu*(Zmean-gridZ)/dZ + elem_r; % this one is a vector
%         
%         % Fill Ass collection
%         for ik=1:intK
%             for iz=1:intZ
%                 A11          = spdiags(elem_b(:,1,ik,iz),0,inta,inta) + spdiags(elem_a(2:inta,1,ik,iz),-1,inta,inta) + spdiags([0;elem_e(1:inta-1,1,ik,iz)],1,inta,inta);
%                 A12          = spdiags(elem_b(:,2,ik,iz),0,inta,inta) + spdiags(elem_a(2:inta,2,ik,iz),-1,inta,inta) + spdiags([0;elem_e(1:inta-1,2,ik,iz)],1,inta,inta);
%                 Ass{ik,iz}    = [A11,sparse(inta, inta); sparse(inta, inta),A12] + Aswitch;
%             end
%         end
%         
%         % Fill A2 collection
%         jumper =inta*intx;
%         for iz = 1:intZ
%             for ik = 1:intK
%                 A2{iz}((ik - 1) * jumper + 1:ik * jumper,(ik - 1) * jumper + 1:ik * jumper) = Ass{ik,iz};
%             end
%         end
%         for iz = 1:intZ
%             for ik = 1:intK - 1
%                 A2{iz}((ik - 1)*jumper+1:ik*jumper, ik*jumper+1:(ik + 1)*jumper) = (Kdot(ik,iz)/dK) * speye(inta*intx,inta*intx);
%             end
%             A2{iz}((intK - 1)*jumper+1:intK*jumper, (intK - 1)*jumper + 1:intK*jumper)=A2{iz}((intK-1)*jumper+1:intK*jumper,(intK-1)*jumper+1:intK*jumper)+(Kdot(intK,iz)/dK)*speye(inta*intx,inta*intx);
%         end
%         
%        % Fill A3
%         jumper = inta*intx*intK;
%         for iz=1:intZ
%             WW((iz-1)*jumper+1:iz*jumper,(iz-1)*jumper+1:iz*jumper) = A2{iz};
%         end
%         for iz=1:intZ-1
%             WW(iz*jumper+1:(iz+1)*jumper,(iz-1)*jumper+1:iz*jumper) = elem_r*speye(inta*intx*intK);
%         end
%         for iz=1:intZ-1
%             WW((iz-1)*jumper+1:iz*jumper,iz*jumper+1:(iz+1)*jumper) = elem_x(iz)*speye(inta*intx*intK);
%         end
%         WW((intZ-1)*jumper+1:intZ*jumper,(intZ-1)*jumper+1:intZ*jumper) = WW((intZ-1)*jumper+1:intZ*jumper,(intZ-1)*jumper+1:intZ*jumper) + elem_x(intZ)*speye(inta*intx*intK, inta*intx*intK);
%         WW(1:jumper,1:jumper) = WW(1:jumper,1:jumper) + elem_r * speye(inta*intx*intK);
% 
%         if isreal(WW) ~= 1
%             disp('')
%             disp('Matrix is not real number')
%             disp('')
%             break
%         end
        
        % Construct sparse matrix for individual wealth
        ElemX = -min(sb,0)/da;
        ElemY = -max(sf,0)/da + min(sb,0)/da;
        ElemZ = max(sf,0)/da;

        % Sparse matrix for indivisdual wealth
        Ass = cell(intK,intZ);
        for iz = 1:intZ
            for ik = 1:intK
                AA1 = spdiags(ElemY(:,1,ik,iz),0,inta,inta) + spdiags(ElemX(2:inta,1,ik,iz),-1,inta,inta) + spdiags([0;ElemZ(1:inta-1,1,ik,iz)],1,inta,inta);
                AA2 = spdiags(ElemY(:,2,ik,iz),0,inta,inta) + spdiags(ElemX(2:inta,2,ik,iz),-1,inta,inta) + spdiags([0;ElemZ(1:inta-1,2,ik,iz)],1,inta,inta);
                AA = [AA1, sparse(inta,inta); sparse(inta,inta), AA2] + Aswitch;

                Ass{ik,iz} = AA + Zcenter{ik,iz} + Kcenter{ik,iz};
            end
        end

        AAss = reshape(Ass, intK*intZ,1);

        % Sparse matrix for calculating value function
        W = cell(intZ*intK,intZ*intK);
        for iw = 1:intZ*intK
            for jw = 1:intZ*intK
                W{iw,jw} = sparse(inta*intx,inta*intx);
            end
        end

        for iw = 1:intZ*intK
            if iw == 1 

                W{iw,iw} = AAss{iw,1};
                W{iw, iw + 1} = KKup{iw,1};
                W{iw, iw + intK} = ZZup{iw,1};

            elseif iw == intZ*intK

                W{iw,iw} = AAss{iw,1};
                W{iw, iw - 1} = KKdown{iw,1};
                W{iw, iw - intK} = ZZdown{iw,1};

            elseif (intK + 1 > iw) && (iw > 1)

                W{iw,iw} = AAss{iw,1};
                W{iw, iw + 1} = KKup{iw,1};
                W{iw, iw - 1} = KKdown{iw,1};
                W{iw, iw + intK} = ZZup{iw,1};

            elseif (intZ*intK > iw) && (iw > intK*(intZ - 1))

                W{iw,iw} = AAss{iw,1};
                W{iw, iw + 1} = KKup{iw,1};
                W{iw, iw - 1} = KKdown{iw,1};
                W{iw, iw - intK} = ZZdown{iw,1};

            else

                W{iw,iw} = AAss{iw,1};
                W{iw, iw + 1} = KKup{iw,1};
                W{iw, iw - 1} = KKdown{iw,1};
                W{iw, iw + intK} = ZZup{iw,1};
                W{iw, iw - intK} = ZZdown{iw,1};

            end
        end
            
        % Translate for Matrix W and debug for WW
        WW = cell2mat(W);
%         spy(WW,'b')
        
        % Calculate Coefficient
        BB = (rho + 1/Delta) * speye(inta*intx*intK*intZ) - WW;
        u_stacked = uss(:); V_stacked = Vss(:);

        bb = u_stacked + V_stacked/Delta;

        Vnew_stacked = BB\bb; %SOLVE SYSTEM OF EQUATIONS
        
        vss = reshape(Vnew_stacked,size(Vss));
        Vchange = Vnew_stacked - V_stacked;
%         Vss = vss;
        
        dist(n) = max(abs(Vchange));
        
        disp('Value function calculating...')
        disp(dist(n));
        
        if dist(n) < crit
            disp('Value Function Converged, Iteration = ')
            disp(n)
            
            % check that the algorithm has converged correctly

            HJB_check = zeros(inta, intx, intK, intZ);
            if gamma == 1
                HJB_check(:,1,:,:) = -rho*vss(:,1,:,:) + log(cs(:,1,:,:)) + (w(:,1,:,:)* (1 - tau).*quadx(:,1,:,:) + w(:,1,:,:)* com.*(1 - quadx(:,1,:,:)) + r(:,1,:,:).*quada(:,1,:,:)-cs(:,1,:,:)).*Va_Upwind(:,1,:,:) + la(1)*(vss(:,2,:,:)-vss(:,1,:,:)) + PLM(:,1,:,:).*dVdK(:,1,:,:) + mu*(Zmean-quadZ(:,1,:,:)).*dVdZ(:,1,:,:) + ((sigma^2)/2)*dVddZ(:,1,:,:);
                HJB_check(:,2,:,:) = -rho*vss(:,2,:,:) + log(cs(:,2,:,:)) + (w(:,2,:,:)* (1 - tau).*quadx(:,2,:,:) + w(:,2,:,:)* com.*(1 - quadx(:,2,:,:)) + r(:,2,:,:).*quada(:,2,:,:)-cs(:,2,:,:)).*Va_Upwind(:,2,:,:) + la(2)*(vss(:,1,:,:)-vss(:,2,:,:)) + PLM(:,2,:,:).*dVdK(:,2,:,:) + mu*(Zmean-quadZ(:,2,:,:)).*dVdZ(:,2,:,:) + ((sigma^2)/2)*dVddZ(:,2,:,:);
            else
                HJB_check(:,1,:,:) = -rho*vss(:,1,:,:) + (cs(:,1,:,:).^(1-gamma))/(1-gamma) + (w(:,1,:,:)* (1 - tau).*quadx(:,1,:,:) + w(:,1,:,:)* com.*(1 - quadx(:,1,:,:)) + r(:,1,:,:).*quada(:,1,:,:)-cs(:,1,:,:)).*Va_Upwind(:,1,:,:) + la(1)*(vss(:,2,:,:)-vss(:,1,:,:)) + PLM(:,1,:,:).*dVdK(:,1,:,:) + mu*(Zmean-quadZ(:,1,:,:)).*dVdZ(:,1,:,:) + ((sigma^2)/2)*dVddZ(:,1,:,:);
                HJB_check(:,2,:,:) = -rho*vss(:,2,:,:) + (cs(:,2,:,:).^(1-gamma))/(1-gamma) + (w(:,2,:,:)* (1 - tau).*quadx(:,2,:,:) + w(:,2,:,:)* com.*(1 - quadx(:,2,:,:)) + r(:,2,:,:).*quada(:,2,:,:)-cs(:,2,:,:)).*Va_Upwind(:,2,:,:) + la(2)*(vss(:,1,:,:)-vss(:,2,:,:)) + PLM(:,2,:,:).*dVdK(:,2,:,:) + mu*(Zmean-quadZ(:,2,:,:)).*dVdZ(:,2,:,:) + ((sigma^2)/2)*dVddZ(:,2,:,:);
            end
%             HJB_check(:,1,:,:) = -rho*vss(:,1,:,:) + (cs(:,1,:,:).^(1-gamma)-1)/(1-gamma) + (w(:,1,:,:).*quadx(:,1,:,:)+r(:,1,:,:).*quada(:,1,:,:)-cs(:,1,:,:)).*Va_Upwind(:,1,:,:) + la(1)*(vss(:,2,:,:)-vss(:,1,:,:)) + PLM(:,1,:,:).*dVdK(:,1,:,:) + mu*(Zmean-quadZ(:,1,:,:)).*dVdZ(:,1,:,:) + ((sigma^2)/2)*dVddZ(:,1,:,:);
% 
%             HJB_check(:,2,:,:) = -rho*vss(:,2,:,:) + (cs(:,2,:,:).^(1-gamma)-1)/(1-gamma) + (w(:,2,:,:).*quadx(:,2,:,:)+r(:,2,:,:).*quada(:,2,:,:)-cs(:,2,:,:)).*Va_Upwind(:,2,:,:) + la(2)*(vss(:,1,:,:)-vss(:,2,:,:)) + PLM(:,2,:,:).*dVdK(:,2,:,:) + mu*(Zmean-quadZ(:,2,:,:)).*dVdZ(:,2,:,:) + ((sigma^2)/2)*dVddZ(:,2,:,:);
            
            HJB_check_stacked = zeros(inta*intx*intK*intZ,1);
            for iZ=2:intZ-1    % reshape could be faster, but this is clearer
                for iK=2:intK-1
                    for ix=1:intx
                        HJB_check_stacked(2+(ix-1)*inta+(iK-1)*inta*intx+(iZ-1)*inta*intx*intK:(inta-1)+(ix-1)*inta+(iK-1)*inta*intx+(iZ-1)*inta*intx*intK,1)=HJB_check(2:inta-1,ix,iK,iZ);
                    end
                end
            end

            disp('HJB_check:')
            disp(max(max(abs(HJB_check_stacked))))
            break
        end
    end
    
    % Calculate Consumption function and Policy function
    % Forward Difference : Indivudual wealth
    Vsaf(1:inta-1,:,:,:) = (Vss(2:inta,:,:,:)-Vss(1:inta-1,:,:,:))/da;
    Vsaf(inta,:,:,:) = (w(inta,:,:,:) * (1 - tau).*quadx(inta,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(inta,:,:,:))+ r(inta,:,:,:).*amax).^(-gamma);
%         Vsaf(inta,:,:,:) = (w(inta,:,:,:).*quadx(inta,:,:,:) + r(inta,:,:,:).*amax).^(-gamma);

    % Backward Difference : Indivisual wealth
    Vsab(2:inta,:,:,:) = (Vss(2:inta,:,:,:)-Vss(1:inta-1,:,:,:))/da; %w(inta,:,:,:) * (1 - tau).*quadx(inta,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(inta,:,:,:))
    Vsab(1,:,:,:) = (w(1,:,:,:) * (1 - tau).*quadx(1,:,:,:) + w(inta,:,:,:) * com.*(1 - quadx(1,:,:,:)) + r(1,:,:,:).*amin).^(-gamma);
%         Vsab(1,:,:,:) = (w(1,:,:,:).*quadx(1,:,:,:) + r(1,:,:,:).*amin).^(-gamma);

    % Consumption and savings with forward difference
    cf = (Vsaf).^(-1/gamma); 
    sf = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cf;
%         sf = w.*quadx + r.*quada - cf;

    % Consumption and savings with backward difference
    cb = (Vsab).^(-1/gamma);
    sb = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cb;
%         sb = w.*quadx + r.*quada - cb;

    % Consumption and derivative of value function at steady state
    c0 = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada;
    Va0 = c0.^(-gamma);
%         c0 = w.*quadx + r.*quada;

    % dV_upwind makes a choice of forward or backward differences based on
    % The sign of the drift for individual wealth
    Iaf = (sf > 0);             % Positive drift --> Forward difference
    Iab = (sb < 0);           % Negative drift --> Backward difference
    Ia0 = (1 - Iaf - Iab); % At steady state

    Va_Upwind = Vsaf.*Iaf + Vsab.*Iab + Va0.*Ia0; %important to include third term
    cs = (Va_Upwind).^(-1/gamma); ps = w * (1 - tau).*quadx + w * com.*(1 - quadx) + r.*quada - cs; 
    if gamma == 1
        uss = log(cs);
    else
        uss = (cs.^(1-gamma) - 1)/(1-gamma);
    end
end