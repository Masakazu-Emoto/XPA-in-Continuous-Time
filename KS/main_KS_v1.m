%% main_KS.m : This code solves the Krusell and Smith model in continuous time by KS Algorithm
%
% Masakazu Emoto and Takeki Sunakawa (2021)
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time"
%
% Reference : Fernandez-Villaverde, J., S. Hurtado and G. Nuno (2019, FVHN hereafter)
% "Solving the Krusell-Smith (1998) model"
%
% The original code is downloaded from https://github.com/jesusfv/financial-frictions/tree/master/KS_LR
%
% Author : Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Revised by Takeki Sunakawa 2021/01/05
% E-mail address : masakazu.emoto@gmail.com
%
% Uses : steadystate.m, inner_v1.m, fokker_planck_v1,m, simulate_v1.m

%% Summary of the algorithm
% NOTE: Steps 0,1 and 2-1 are common between KS and XPA algorithms
% Step 0 : Set parameters
% Step 1 : Solve for the deterministic steady state
% Step 2 : KS algorithm
% Step 2-1 : Inner Loop, Calculate the policy function by taking the forecasting rule (perceived law of motion) as given 
% Step 2-2 : Outer Loop (1), Simulate the path of aggregate capital
% Step 2-3 : Outer Loop (2), Solve for the forecasting rule by linear
% regression of simulated data
% (Step 3 : Solve for the stochastic steady state)
% Step 4 : Simulate the model and calculate the Den Haan Error
% Step 5 : Plot graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; format long;
clc; tic;

%% NOTE: This code is based on the ones written by FVHN. However, we extend their original code in the following two dimensions:
UpwindKZ = 1; % (1) We use the upwind scheme not only individual wealth, a, but also K and Z.
KFEnoKZ  = 1; % (2) We exclude the direct effect of aggregate variables K and Z on the matrix
% A_lm in their note when we solve the KF equation (there is the indirect effect through r and w).

%% -------------------------------------------------- %
%% Step 0 : Set Parameters
%% -------------------------------------------------- %
parameters;

%% -------------------------------------------------- %
%% Step 1 : Solve for the deterministic steady state
%% -------------------------------------------------- %
% Calculate deterministic steady state
% tic;
disp('Calcutating the deterinistic steady state')
[rds, wds, Kds, Ads, uds, cds, pds, ids, Vds, gds] = steadystate();
% toc;

% disp('Deterministic steady state')
disp(Kds) 
% disp(sum(gds'*grida*da))

% Set the grid around the deterministic steady state of K
if sigma >= 0.03
    Kmax = 1.15*Kds; Kmin = 0.85*Kds; %intK = 3;
else
    Kmax = 1.05*Kds; Kmin = 0.95*Kds; %intK = 3;
end
gridK = linspace(Kmin,Kmax,intK)'; dK = (Kmax - Kmin)/(intK - 1);

%% -------------------------------------------------- %
%% Step 2 : KS algorithm
%% -------------------------------------------------- %
% Resize grid
global quada quadx quadK quadZ
quada = zeros(inta,intx,intK,intZ);
quadx = zeros(inta,intx,intK,intZ);
quadK = zeros(inta,intx,intK,intZ);
quadZ = zeros(inta,intx,intK,intZ);

for ix=1:intx
    for ik=1:intK
        for iz=1:intZ
            quada(:,ix,ik,iz)=grida;
        end
    end
end

for ia=1:inta
    for ik=1:intK
        for iz=1:intZ
            quadx(ia,:,ik,iz)=x;
        end
    end
end

for ia=1:inta
    for ix=1:intx
        for iz=1:intZ
            quadK(ia,ix,:,iz)=gridK;
        end
    end
end

for ia=1:inta
    for ix=1:intx
        for ik=1:intK
            quadZ(ia,ix,ik,:)=gridZ;
        end
    end
end

% Initial guess for the forecasting rule: the same as in FVHN
Kdot = zeros(intK,intZ); 

% Initial guess of vss, r, and w for the inner loop (used in inner.m)
r = alpha*quadK.^(alpha - 1).*LAve^(1 - alpha).*exp(quadZ) - delta;
w = (1 - alpha)*quadK.^(alpha).*LAve^(-alpha).*exp(quadZ);
if gamma == 1
    vss = log((w.*(1 - tau).*quadx + w.*com.*(1 - quadx) + r.*quada))/rho;
else
    vss = (w.*(1 - tau).*quadx + w.*com.*(1 - quadx) + r.*quada).^(1-gamma)/(1-gamma)/rho;
end

% Initial distribution from the deterministic steady state
muini = gds;

% Inner loop and outer loop
disp('Calculating the inner and outer loops by KS algorithm')
disp(' ')
for iteration=1:maxitK % Outer loop
%while (epsilon > epsmin)
    %% -------------------------------------------------- %
    %% Step 2-1 : Inner Loop, Calculate the policy function by taking the forecasting rule (perceived law of motion) as given 
    %% -------------------------------------------------- %
%     tic;
%    [A1, A1tilde, A3, vss, cs, ps] = inner(Kdot, vss, iteration, r, w);
    [A1, A1tilde, A3, vss, cs, ps] = inner_v1(Kdot, vss, r, w, UpwindKZ);
    disp('  Finished solving the HJB Equation')
%     toc;
    
    %% -------------------------------------------------- %
    %% Step 2-2 : Outer Loop (1), Simulate the path of aggregate capital
    %% -------------------------------------------------- %
%     tic;
    if (KFEnoKZ)
        [Ksim, Zsim, Kdown, Kup, Zdown, Zup] = fokker_planck_v1(Zshocks, muini, A1tilde);
    else
        [Ksim, Zsim, Kdown, Kup, Zdown, Zup] = fokker_planck_v1(Zshocks, muini, A1);
    end
    disp('  Finished simulating aggregate capital')
%     toc;
    
    %% -------------------------------------------------- %
    %% Step 2-3 : Outer Loop (2), Solve for the forecasting rule by linear
    %% regression of simulated data
    %% -------------------------------------------------- %
    % linear regression
    Y = (Ksim(Dtime+1:Stime) - Ksim(Dtime:Stime-1))/dT; % dK(t): Growth rate of K at time t
    X0 = ones(Stime-Dtime,1);                           % Constant
    X1 = Ksim(Dtime:Stime-1);                           % K(t-1)
    X2 = Zsim(Dtime:Stime-1);                           % Z(t-1)

    X1 = log(X1);
    %% NOTE: FVHN use cross terms, whereas we don't
    %X3 = X1.*X2;   % log(K(t-1))*Z(t-1)
    X = [X0 X1 X2];
    %X = [X0 X1 X2 X3];
    
%    B = (X'*X)^-1*X'*Y;
    B = (X'*X)\(X'*Y);
    Y_LR = X * B;
    Y_Stad = (mean(Y - Y_LR).^2).^0.5;
    Y_R2 = 1 - (sum((Y - Y_LR).^2))/(sum((Y - mean(Y)).^2));
    
    % Use the coefficients to calculate the forecasting rule on HJB grid
    X1mm = squeeze(quadK(1,1,:,:)) ;
    X2mm = squeeze(quadZ(1,1,:,:)) ;

    X1m = reshape(X1mm,[intK*intZ,1]);
    X2m = reshape(X2mm,[intK*intZ,1]);
    X0m = ones(size(X1m));
    
    X1m=log(X1m); 
    %X3m = X1m.*X2m;
    X_LR = [X0m X1m X2m];
    %X_LR = [X0m X1m X2m X3m];
    
    Kdotnew = X_LR*B;
    Kdotnew = reshape(Kdotnew, size(Kdot)); % Same value ! ????????
    
    if ~isreal(Kdotnew)
        disp('')
        disp('  The matrix for aggregate dynamics is not real')
        disp('')
        break
    end
    
    epsilon = max(max(abs(Kdot - Kdotnew)));
    
    if epsilon > critK
%    if epsilon > epsmin
%        disp('Calculating the law of motion ...')
%        disp([Y_R2, epsilon])
%        disp([iteration, epsilon, Y_R2])
        fprintf("  iter = %4d, diff = %5.6f, R2 = %5.6f\n",iteration, epsilon, Y_R2)
    else
        break;
    end
    
    % Update law of motion
    % TS: Why do we do this?
    Kdot = (1 - relax_dot) * Kdot + relax_dot * Kdotnew;
    relax_dot = relax_dot * relax1 + relax2;

end

toc;

%% -------------------------------------------------- %
%% Step 3 : Solve for the stochastic steady state
%% -------------------------------------------------- %
% disp('Calculating stochastic steady state')
% muini = gds; Zss = ones(N,1).*Zmean;
% ssTFP = exp(Zss);
% [ss_mu, Kss, KKss] = stochastic_steady(simZ, muini, Bss, Kdotnew);
% toc;

%% -------------------------------------------------- %
%% Step 4 : Simulate the model and calculate the Den Haan Error
%% -------------------------------------------------- %
disp('Simulating the model and Calculating Den Haan Error')
% We refer to Ahn et al to calculate the Den Haan error
N = 10000; 
muini = gds; 
Zsim = zeros(N,1); 
rng(100);  
shock = randn(N,1); 
shock(1,1) = 0; % why ???
mmu = -1 + mu;
for time = 1:N-1
    if time == 1
        Zsim(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT));
    else
        Zsim(time+1) = (1 - mmu * dT)^(-1) * (Zsim(time) + sigma * shock(time) * sqrt(dT));
    end
end
simTFP = exp(Zsim);
if (KFEnoKZ)
    [sim_mu, KS_K, KS_KK] =  simulate_v1(Zsim, muini, A1tilde, Kdotnew);
else
    [sim_mu, KS_K, KS_KK] =  simulate_v1(Zsim, muini, A1, Kdotnew);
end

% KS_KK is the sequence of simulated results using the forecasting rule only
% KS_K is the sequence of simulated results using the forecasting rule and the HJB equation
Drop = 1000; 
DH_Error = 100.0 * max(abs(log(KS_KK(Drop+1:end)) - log(KS_K(Drop+1:end))));
DH_Mean = 100.0 * sum(abs(log(KS_KK(Drop+1:end)) - log(KS_K(Drop+1:end))))/(N - Drop);

disp('MAX Den Haan Error')
disp(DH_Error)
disp('MEAN Den Haan Error')
disp(DH_Mean)

toc;

% disp('Stochastic steady state by KS Algorithm')
% disp(Kss(end)) 

%% -------------------------------------------------- %
%% Step 5 : Plot graphs
%% -------------------------------------------------- %
figure(1)
surf(gridZ,gridK,Kdot);
title('The perceived law of motion, PLM : KS', 'interpreter','latex','FontSize',14);
xlabel('shock ($Z$)', 'interpreter','latex','FontSize',14);
ylabel('capital ($K$)', 'interpreter','latex','FontSize',14);
xlim([Zmin Zmax]); ylim([Kmin Kmax]);

intKK = 16; 
intZZ = 41; 
gridKK = linspace(Kmin,Kmax,intKK)';
gridZZ = linspace(Zmin, Zmax,intZZ)';
LOM = interp1(gridK,Kdot,gridKK,'spline');
LOM = interp1(gridZ,LOM',gridZZ,'spline');

figure(2)
surf(gridZZ,gridKK,LOM');
title('Forecasting rule : KS : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Shock : $Z$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
zlabel('$\dot{K}$', 'interpreter','latex','FontSize',10);
xlim([Zmin Zmax]); ylim([Kmin Kmax]);

% Ploting Sparse Matrix for value function
figure(3)
spy(A3,'b');    

% Plot Transition Dynamics from DSS to SSS
% figure(4)
% plot(Kss,'b-','LineWidth',1); grid;
% hold on
% plot(Kss(1,1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1); grid;
% hold on
% plot(max(size(Kss)),Kss(end,1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b');
% title('Transition Dynamics the DSS to the SSS : $\sigma$ = 0.07', 'interpreter','latex','FontSize',10);
% xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
% ylabel('Capital : $K$', 'interpreter','latex','FontSize',10); grid;

% Plotting Simulaion Path for aggregate capital
figure(5)
plot(KS_K,'b-','LineWidth',1);
hold on
plot(KS_KK,'r-','LineWidth',1); 
title('Simulaton Path : KS : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('$K^*_{t}$', '$\tilde{K}_{t}$','Location','northwest','interpreter','latex'); grid;
toc;

% Save the Results
% save CT_KS.mat