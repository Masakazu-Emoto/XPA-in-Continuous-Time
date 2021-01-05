%% main_XPA.m : This code solves the Krusell and Smith model in continuous time by XPA Algorithm
%
% Masakazu Emoto and Takeki Sunakawa (2021)
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time"
%
% Reference : Fernandez-Villaverde, J., S. Hurtado and G. Nuno (2019, FVHN hereafter)
% "Solving the Krusell-Smith (1998) model" "Financial Frictions and the Wealth Distribution"
% Sunakawa, T. (2020, Computational Economics)
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Macro Models"
%
% The original code is downloaded from 
% https://github.com/jesusfv/financial-frictions/tree/master/KS_LR
% https://github.com/tkksnk/Xpa
%
% Author : Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Revised by Takeki Sunakawa 2021/01/05
% E-mail address : masakazu.emoto@gmail.com
%
% Uses : steadystate.m, inner_v1.m, bias.m, outer.m, simulate_v1.m

%% Summary of the algorithm
% NOTE: Steps 0,1 and 2-1 are common between KS and XPA algorithms
% Step 0 : Set Parameters
% Step 1 : Solve for the deterministic steady state
% Step 2 : XPA algorithm
% Step 2-0 : Solve for bias correction terms
% Step 2-1 : Inner Loop, Calculate the policy function by taking the forecasting rule (perceived law of motion) as given 
% Step 2-2 : Outer Loop, Calculate the forecasting rule by taking the policy function as given
% (Step 3 : Solve for the stochastic steady state)
% Step 4 : Simulate the model and calculate the Den Haan Error
% Step 5 : Plot graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; format long;
clc; tic;
addpath ../common

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
% Calculate deterministic steady state (Deterministic steady state is the steady state without aggregate uncertainty)
% tic;
disp('Calcutating deterinistic steady state')
[rds, wds, Kds, Ads, uds, cds, pds, ids, Vds, gds, X, Y, Z] = steadystate();
% toc;

disp(Kds) 
% disp(sum(gds'*grida*da))

% Set the grid around the deterministic steady state of K
if sigma >= 0.03
    Kmax = 1.15*Kds; Kmin = 0.85*Kds;
else
    Kmax = 1.05*Kds; Kmin = 0.95*Kds;
end
gridK = linspace(Kmin,Kmax,intK)'; dK = (Kmax - Kmin)/(intK - 1);

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

%% -------------------------------------------------- %
%% Step 2 : XPA algorithm
%% -------------------------------------------------- %

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

% Inner loop and outer loop
disp('Calcutating inner loop and outer loop by XPA algorithm')
disp(' ')
%% -------------------------------------------------- %
%% Step 2-0 : Solve for bias correction terms
%% -------------------------------------------------- %
[muxz, psix, zeta] = bias(gds, pds);

for iteration=1:maxitK % Outer loop
    %% -------------------------------------------------- %
    %% Step 2-1 : Inner Loop, Calculate the policy function by taking the forecasting rule (perceived law of motion) as given 
    %% -------------------------------------------------- %
%     tic;
%    [A1, A1tilde, A3, vss, cs, ps] = inner(Kdot, vss, iteration, r, w);
    [A1, A1tilde, A3, vss, cs, ps] = inner_v1(Kdot, vss, r, w, UpwindKZ);
    disp('  Finished solving the HJB Equation')
%     toc;
    
    %% -------------------------------------------------- %
    %% Step 2-2 : Outer Loop, Calculate the forecasting rule by taking the policy function as given
    %% -------------------------------------------------- %
    Kdotnew = outer(ps, muxz, psix, zeta);
    disp('  Finished calculating the forecasting rule')
    
    if ~isreal(Kdotnew)
        disp('')
        disp('  The matrix for aggregate dynamics is not real')
        disp('')
        break
    end
    
    epsilon = max(max(abs(Kdot - Kdotnew)));
    
    if epsilon > critK
        fprintf("  iter = %4d, diff = %5.6f\n",iteration, epsilon)
    else
        break;
    end
    
    % Update the forecasting rule
    Kdot = relax * Kdot + (1 - relax) * Kdotnew;
    
end

toc;

%% -------------------------------------------------- %
%% Step 3 : Solve for the stochastic steady state
%% -------------------------------------------------- %
% disp('Calculating stochastic steady state')
% muini = gds; Zss = ones(N,1).*Zmean;
% ssTFP = exp(Zss);
% [ss_mu, Kss, KKss] = stochastic_steady(Zss, muini, Bss, Kdot);
% toc;

% disp('Stochastic steady state by XPA Algorithm')
% disp(Kss(end))

%% -------------------------------------------------- %
%% Step 4 : Simulate the model and Calcuate the Den Haan Error
%% -------------------------------------------------- %
disp('Simulating the model and Calculating Den Haan Error')
% We refer to Ahn et al to calculate the Den Haan error
N = 10000; 
muini = gds; 
simZ = zeros(N,1); 
rng(100); % Simulate Periods 125
shock = randn(N,1); 
%shock(1,1) = 0; 
mmu = -1 + mu; 
for time = 1:N-1
    if time == 1
        simZ(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT));
    else
        simZ(time+1) = (1 - mmu * dT)^(-1) * (simZ(time)  + sigma * shock(time) * sqrt(dT));
    end
end

simTFP = exp(simZ);
if (KFEnoKZ)
    [sim_mu, XPA_K, XPA_KK] =  simulate_v1(simZ, muini, A1tilde, Kdotnew);
else
    [sim_mu, XPA_K, XPA_KK] =  simulate_v1(simZ, muini, A1, Kdotnew);
end

% XPA_KK is the sequence of simulated results using the forecasting rule only
% XPA_K is the sequence of simulated results using the forecasting rule and the HJB equation
Drop = 1000;
DH_Error = 100.0 * max(abs(log(XPA_KK(Drop+1:end)) - log(XPA_K(Drop+1:end))));
DH_Mean = 100.0 * sum(abs(log(XPA_KK(Drop+1:end)) - log(XPA_K(Drop+1:end))))/(N - Drop);

disp('MAX Den Haan Error')
disp(DH_Error)
disp('MEAN Den Haan Error')
disp(DH_Mean)

toc;

%% -------------------------------------------------- %
%% Step 5 : Plot graphs
%% -------------------------------------------------- %
figure(1)
surf(gridZ,gridK,Kdot);
title('The aggregate law of motion, PLM : XPA', 'interpreter','latex','FontSize',10);
xlabel('Shock ($Z$)', 'interpreter','latex','FontSize',10);
ylabel('Capital ($K$)', 'interpreter','latex','FontSize',10);
xlim([Zmin Zmax]); ylim([Kmin Kmax]);

intKK = 16; 
intZZ = 41; 
gridKK = linspace(Kmin,Kmax,intKK)';
gridZZ = linspace(Zmin, Zmax,intZZ)';
LOM = interp1(gridK,Kdot,gridKK,'spline');
LOM = interp1(gridZ,LOM',gridZZ,'spline');

figure(2)
surf(gridZZ,gridKK,LOM');
title('Forecasting rule : XPA : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Shock : $Z$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
zlabel('$\dot{K}$', 'interpreter','latex','FontSize',10);
xlim([Zmin Zmax]); ylim([Kmin Kmax]);

% Ploting Sparse Matrix for value function
figure(3)
spy(A3,'b');    

% Plotting Transition Dynamics DSS to SSS
% figure(4)
% plot(Kss,'b-','LineWidth',1);
% hold on
% plot(Kss(1,1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
% hold on
% plot(max(size(Kss)),Kss(end,1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b');
% title('Transition Dynamics the DSS to the SSS : $\sigma$ = 0.07', 'interpreter','latex','FontSize',10);
% xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
% ylabel('Capital : $K$', 'interpreter','latex','FontSize',10); grid;

% Plotting Simulaion Path for aggregate capital
figure(5)
plot(XPA_K,'b-','LineWidth',1);
hold on
plot(XPA_KK,'r-','LineWidth',1); 
hold on
title('Simulaton Path : XPA : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('$K^*_{t}$', '$\tilde{K}_{t}$','Location','northwest','interpreter','latex'); grid;
xlim([1 N]);

% Save the results
% save CT_XPA.mat