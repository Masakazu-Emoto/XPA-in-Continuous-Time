%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time."
% By Masakazu Emoto and Takeki Sunakawa
% This code solves the Krusell and Smith model in continuous time by XPA Algorithm
% XPA algoithm is cited by Sunakawa (2020 Computational Economics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm 
% Step 0 : Set Parameters
% Step 1 : Solve for deterministic steady state
% Step 2 : Solve for bias correction terms
% Step 3-1 : Inner Loop Taking as given the forecasting rule, calculating policy function
% Step 3-2 : Outer Loop Taking as given policy function, calculating the forecasting rule
% (Step 4 : Solve for stochastic steady state)
% Step 5 : Simulation the model and Calculate the Den Haan Error
% Step 6 : Plot relevant values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; format long;
clc; tic;

% -------------------------------------------------- %
% Step 0 : Set Parameters
% -------------------------------------------------- %
global gamma rho alpha delta la intx x mu sigma com tau LAve 
global maxit maxitK crit critK Delta damp relax

% Preference
gamma = 1;   % Coefficient of relative risk aversion
rho = 0.01;     % Discount rate

% Production
alpha = 0.36;    % Capital share
delta = 0.025;  % Depriciate rate

% Idiosyncratic shock for labor productivity
intx = 2;
x1 = 0; x2 = 1;
x = [x1, x2]; % Labpr producticvity x(1) : unemployed x(2) :employed

% Transition rate for labor productivity 
% lambda1 : unemployment -> employment
% lambda2 : employment -> unemployment
lambda1 = 0.5;                                         
lambda2 = (lambda1 / (x(2) * 0.93 - x(1))*(x(2) - x(2) * 0.93));
la = [lambda1, lambda2];

% Aggregate shock for TFP (Ornstein-uhlenbeck process)
mu = 0.25;       % Mean
sigma = 0.007; % Variance

% Tax system (Ahn et al.(2018))
com = 0.15; 
tau = (com / x(2)) * (la(2) / la(1));

% Average labour supply
LAve = (la(1) * x(2) + la(2) * x(1)) / (la(1) + la(2));

% PARAMETERS : Convergence Criterion is refered by Villaverde et al.(2019)
maxit  = 100;    % maximum number of iterations in the HJB loop
maxitK = 100;   % maximum number of iterations in the K loop
crit = 1e-6;        % criterion HJB loop
critK = 1e-5;     % criterion K loop
Delta = 1000;    % delta in HJB algorithm
damp = 0.001;  % relaxation parameter for Steady state
relax = 0.9;        % relaxation parameter for Law of motion

% -------------------------------------------------- %
% Setting grid
% -------------------------------------------------- %
global inta amin amax grida da aa aaa xx xxx Aswitch

% Indivisual wealth grid
inta = 100; amin = 0; amax = 100;
grida = linspace(amin,amax,inta)';
da = (amax - amin) / (inta - 1);
aa = [grida,grida];
aaa = reshape(aa,2*inta,1);

% Labor productivity grid
xx = ones(inta,1) * x;
xxx = reshape(xx,2*inta,1);

% Idiosyncratic shocks for labor productivity
Aswitch = [-speye(inta) * la(1), speye(inta) * la(1); speye(inta) * la(2), -speye(inta) * la(2)];

% -------------------------------------------------- %
% Step 1 : Solve for deterministic steady state
% -------------------------------------------------- %

% Calculate deterministic steady state (Deterministic steady state is the steady state without aggregate uncertainty)
disp('Calcutating deterinistic steady state')
[rds, wds, Kds, Ads, uds, cds, pds, ids, Vds, gds, X, Y, Z] = steadystate();
toc;

disp('Deterministic steady state')
disp(Kds) 
disp(sum(gds'*grida*da))

% -------------------------------------------------- %
% Step 2 : Solve for bias correction terms
% -------------------------------------------------- %

% Setting aggregate wealth and productivity grid
global Kmax Kmin intK gridK dK
global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ

% COMMENT : What is the optimum number of grids?
if sigma >= 0.03
    Kmax = 1.15*Kds; Kmin = 0.85*Kds; intK = 5;
else
    Kmax = 1.05*Kds; Kmin = 0.95*Kds; intK = 5;
end
gridK = linspace(Kmin,Kmax,intK)'; dK = (Kmax - Kmin)/(intK - 1);

% COMMENT : What is the optimum number of grids?
Zmax = 2*sigma; Zmin = -2*sigma; intZ = 5; Zmean = 0; 
gridZ = linspace(Zmin,Zmax,intZ)'; dZ = (Zmax - Zmin)/(intZ - 1); ddZ = dZ^2;
gridZ((intZ+1)/2,1) = Zmean;

% Aggregate Shock Process
zmu = mu.*(Zmean - gridZ); 
zsigma = sigma^2.*ones(intZ,1);

% Calculate the bias correction terms
[zeta, psix, muxz, mpvec] = bias(gds, pds);

% -------------------------------------------------- %
% Step 3 : XPA algorithm
% -------------------------------------------------- %
disp('Calcutating law of motion by XPA algorithm')

% Initial guess for law of motion for Krusell and Smith (1998)
% Initial guess is the same as Villaverde et al.(2019)
Kdot = zeros(intK,intZ);

disp('Calcutating inner loop and outer loop by XPA algorithm')
epsmin = 1e-6; epsilon = 1e+6; iteration = 1;
global  quada quadx quadK quadZ

% Resize grid
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

% Initial guess for XPA Algorithm
r = alpha*quadK.^(alpha - 1).*LAve^(1 - alpha).*exp(quadZ) - delta;
w = (1 - alpha)*quadK.^(alpha).*LAve^(-alpha).*exp(quadZ);
if gamma == 1
    vss = log((w.*(1 - tau).*quadx + w.*com.*(1 - quadx) + r.*quada))/rho;
else
    vss = (w.*(1 - tau).*quadx + w.*com.*(1 - quadx) + r.*quada).^(1-gamma)/(1-gamma)/rho;
end

while (epsilon > epsmin)
    % -------------------------------------------------- %
    % Step 3-1 : Inner Loop Taking as given the forecasting rule, calculating policy function
    % -------------------------------------------------- %
   [Ass, Bss, WW, vss, cs, ps, zx, zy, zz]   = inner(Kdot, vss, iteration, r, w);
%     [Ass, Bss, WW, vss, cs, ps, zx, zy, zz]   = inner_org(Kdot, vss, iteration, r, w);
    disp('Finished calculation HJB Equation')
    
    % -------------------------------------------------- %
    % Step 3-2 : Outer Loop Taking as given policy function, calculating the forecasting rule
    % -------------------------------------------------- %
    [Kdotnew, xpavec] = outer(ps, muxz, psix, zeta);
    disp('Finished calculation law of motion')
    
    if isreal(Kdotnew) ~= 1
        disp('')
        disp('Aggregate dynamics is not real number')
        disp('')
        break
    end
    
    epsilon = max(max(abs(Kdot - Kdotnew)));
    iteration = iteration + 1;
    
    if epsilon > epsmin
        disp('Law of Motion calculating ...')
        disp(epsilon)
    end
    
    % Update law of motion
    Kdot = relax * Kdot + (1 - relax) * Kdotnew;
end

toc;

% -------------------------------------------------- %
% Step 4 : Solve for stochastic steady state
% -------------------------------------------------- %
global T N vtime dT
T = 500; N = 2000; vtime = linspace(0,T,N); % Simulate Periods 200
dT = T/N;

% disp('Calculating stochastic steady state')
% muini = gds; Zss = ones(N,1).*Zmean;
% ssTFP = exp(Zss);
% [ss_mu, Kss, KKss] = stochastic_steady(Zss, muini, Bss, Kdot);
% toc;

% -------------------------------------------------- %
% Step 5 : Simulate the model and Calcuate the Den Haan Error
% -------------------------------------------------- %
disp('Simulating the model')
N = 10000; muini = gds; simZ = zeros(N,1); rng(100);  % Simulate Periods 125
mmu = -1 + mu; shock = randn(N,1); %shock(1,1) = 0; 
for time = 1:N-1
    if time == 1
        simZ(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT)); % Ahn et al method
    else
        simZ(time+1) = (1 - mmu * dT)^(-1) * (simZ(time)  + sigma * shock(time) * sqrt(dT)); % Ahn et al method
    end
end

% XPA_KK is simulated results using the dynamics of aggregate capital
% XPA_K is simulated results using the dynamics of aggregate capital and the HJB equation
simTFP = exp(simZ);
[sim_mu, XPA_K, XPA_KK] =  simulate(simZ, muini, Bss, Kdot);

disp('Calculating Den Haan Error')
Drop = 1000;
DH_Error = 100.0 * max(abs(log(XPA_KK(Drop+1:end)) - log(XPA_K(Drop+1:end))));
DH_Mean = 100.0 * sum(abs(log(XPA_KK(Drop+1:end)) - log(XPA_K(Drop+1:end))))/(N - Drop);

disp('MAX Den Haan Error')
disp(DH_Error)
disp('MEAN Den Haan Error')
disp(DH_Mean)
toc;

% disp('Stochastic steady state by XPA Algorithm')
% disp(Kss(end))

% -------------------------------------------------- %
% Step 6 : Plot relevant values
% -------------------------------------------------- %

% Plotting law of motion
figure(1)
surf(gridZ,gridK,Kdot);
title('The aggregate law of motion, PLM : XPA', 'interpreter','latex','FontSize',10);
xlabel('Shock ($Z$)', 'interpreter','latex','FontSize',10);
ylabel('Capital ($K$)', 'interpreter','latex','FontSize',10);
xlim([Zmin Zmax]); ylim([Kmin Kmax]);

intKK = 16; gridKK = linspace(Kmin,Kmax,intKK)';
intZZ = 41; gridZZ = linspace(Zmin, Zmax,intZZ)';
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
spy(WW,'b');    

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