%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Krusell-Smith Model (1998) in Continuous Time by XPA Algorithm
% XPA algoithm is cited by Sunakawa (2020 Computational Economics)
% Masakazu Emoto @ Kobe univerisity 2020/07/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm 
% Phase 1 : Calculate deterministic steady state
% Phase 2 : Calculate the bias correction terms
% Phase 3-1 : Inner Loop Taking as given the forecasting rule, calculating policy function
% Phase 3-2 : Outer Loop Taking as given policy function, calculating the forecasting rule
% Phase 4 : Calculate stochastic steady state
% Phase 5 : Simulation and Calculate Den Haan Error
% Phase 6 : Plots the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; format long;
clc; tic;

% -------------------------------------------------- %
% Setting Parameter cited by Villaverde et al (2019 NBER)
% -------------------------------------------------- %
global gamma rho alpha delta la intx x mu sigma com tau LAve 
global maxit maxitK crit critK Delta damp relax

% Preference
gamma = 1; % Coefficient of relative risk aversion
rho = 0.01;     % Discount rate

% Production
alpha = 0.36;    % Capital share
delta = 0.025;  % Depriciate rate

% Idiosyncratic shock for labor productivity
intx = 2;
x1 = 0;
x2 = 1;
x = [x1, x2]; % Labpr producticvity x(1) : unemployed x(2) :employed

% Transition rate for labor productivity 
% lambda1 : unemployment -> employment
% lambda2 : employment -> unemployment
lambda1 = 0.5;                                         
lambda2 = (lambda1 / (x(2) * 0.93 - x(1))*(x(2) - x(2) * 0.93));
la = [lambda1, lambda2];

% Aggregate shock for TFP (Ornstein-uhlenbeck process)
mu = 0.25;       % Mean
sigma = 0.05; % Variance 

% Tax system (Ahn et al.(2018))
com = 0.15; 
tau = (com / x(2)) * (la(2) / la(1));

% Average labour supply (Villaverde et al.(2019))
LAve = (la(1) * x(2) + la(2) * x(1)) / (la(1) + la(2));

% PARAMETERS : Convergence Criterion
maxit  = 100;    % maximum number of iterations in the HJB loop
maxitK = 100;   % maximum number of iterations in the K loop
crit = 1e-6;   % criterion HJB loop
critK = 1e-5;    % criterion K loop
Delta = 1000;   % delta in HJB algorithm
damp = 0.001;     % relaxation parameter for Steady state
relax = 0.7;       % relaxation parameter for Law of motion

% -------------------------------------------------- %
% Setting grid cited by Villaverde et al (2019 NBER)
% -------------------------------------------------- %
% Indivisual wealth grid
global inta amin amax grida da aa aaa xx xxx Aswitch

inta = 100;
amin = 0; amax = 100;
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
% Phase 1 : Calculate determinisitec steady state
% -------------------------------------------------- %

% Calculate deterministic steady state
disp('Calcutating deterinistic steady state')
[rds, wds, Kds, Ads, uds, cds, pds, ids, Vds, gds, X, Y, Z] = steadystate();
toc;

disp('Deterministic steady state')
disp(Kds) 
disp(sum(gds'*grida*da))

% Calculate gini coefficient and share of wealth at deterministic steady state
disp('Calcutating gini coefficient at deterinistic steady state')

% -------------------------------------------------- %
% Phase 2 : Calculate the bias correction terms
% -------------------------------------------------- %

% Setting aggregate wealth and productivity grid
global Kmax Kmin intK gridK dK
global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ

% COMMENT : What is the optimum number of grids?
if sigma >= 0.03
    Kmax = 1.15*Kds; Kmin = 0.85*Kds; intK = 3;
else
    Kmax = 1.05*Kds; Kmin = 0.95*Kds; intK = 3;
end
gridK = linspace(Kmin,Kmax,intK)'; dK = (Kmax - Kmin)/(intK - 1);

% COMMENT : What is the optimum number of grids?
Zmax = 2*sigma; Zmin = -2*sigma; intZ = 3; Zmean = 0; 
gridZ = linspace(Zmin,Zmax,intZ)'; dZ = (Zmax - Zmin)/(intZ - 1); ddZ = dZ^2;
gridZ((intZ+1)/2,1) = Zmean;

% Aggregate Shock Process
zmu = mu.*(Zmean - gridZ); 
zsigma = sigma^2.*ones(intZ,1);

% Calculate the bias correction terms
[zeta, psix, muxz, mpvec] = bias(gds, pds);

% -------------------------------------------------- %
% Phase 3 : XPA algorithm
% -------------------------------------------------- %
disp('Calcutating law of motion by XPA algorithm')

% Initial guess for law of motion for Krusell and Smith (1998)
% Initial guess is the same as Villaverde et al.(2019)
Kdot = zeros(intK,intZ);

% Inner loop and outer loop
disp('Calcutating inner loop and outer loop by XPA algorithm')
epsmin = 1e-6; epsilon = 1e+6; iteration = 1;
global  quada quadx quadK quadZ

for ix=1:intx    % repmat would be faster, but this is clearer
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
    % Phase 3-1 : Calculate policy function : Inner Loop
    % -------------------------------------------------- %
    [Ass, WW, vss, cs, ps]   = inner(Kdot, vss, iteration, r, w);
    disp('Finished calculation HJB Equation')
    
    % -------------------------------------------------- %
    % Phase 3-2 : Calculate law of motion : Outer Loop
    % -------------------------------------------------- %
    [Kdotnew, xpavec] = outer(ps, muxz, psix, zeta);
    disp('Finished calculation Forecasting Rule')
    
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
% Phase 4 : Calculate stochastic steady state and Den Haan Error
% -------------------------------------------------- %
global T N vtime dT
T = 500; N = 2000; vtime = linspace(0,T,N); % Simulate Periods 200

dT = T/N;

disp('Calculating stochastic steady state')
muini = gds; simZ = ones(N,1).*Zmean;
simTFP = exp(simZ);
[simvalue, sim_mu, Kss, KKss] = stochastic_steady(simZ, muini, Ass, Kdot);
toc;

% -------------------------------------------------- %
% Phase 5 : Simulation and Calculate Den Haan Error
% -------------------------------------------------- %
disp('Simulating the model')
N = 10000; muini = gds; irfZ = zeros(N,1); rng(100);  % Simulate Periods 125
mmu = -1 + mu; shock = randn(N,1); %shock(1,1) = 0; 
for time = 1:N-1
    if time == 1
        irfZ(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT)); % Ahn et al method
    else
        irfZ(time+1) = (1 - mmu * dT)^(-1) * (irfZ(time)  + sigma * shock(time) * sqrt(dT)); % Ahn et al method
    end
end

% irf_KK : irf_KK is simulated results using the dynamics of aggregate capital
% irf_K : irf_K is simulated results using the dynamics of aggregate capital and the HJB equation
irfTFP = exp(irfZ);
[testvalue, test_mu, XPA_K, XPA_KK] =  simulate(irfZ, muini, Ass, Kdot);

disp('Calculating Den Haan Error')
Drop = 1000;
DH_Error = 100.0 * max(abs(log(XPA_KK(Drop+1:end)) - log(XPA_K(Drop+1:end))));
DH_Mean = 100.0 * sum(abs(log(XPA_KK(Drop+1:end)) - log(XPA_K(Drop+1:end))))/(N - Drop);
DDH_Error = max(abs(XPA_KK(Drop+1:end)-XPA_K(Drop+1:end))./abs(XPA_K(Drop+1:end)));

disp('MAX Den Haan Error')
disp(DH_Error)
disp('MEAN Den Haan Error')
disp(DH_Mean)
toc;

disp('Stochastic steady state by XPA Algorithm')
disp(Kss(end))

% -------------------------------------------------- %
% Phase 6 : Plots the results
% -------------------------------------------------- %
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

% Plot Sparse Matrix
figure(3)
spy(WW,'b');    

% Plot Transition Dynamics DSS to SSS
figure(4)
plot(Kss,'b-','LineWidth',1);
hold on
plot(Kss(1,1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
hold on
plot(max(size(Kss)),Kss(end,1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b');
title('Transition Dynamics the DSS to the SSS : $\sigma$ = 0.07', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10); grid;

% Plot Simulaion Path
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

figure(6)
plot(irfZ,'b-','LineWidth',1.5);
title('Simulaton Path : CT-XPA : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('TFP : $Z$', 'interpreter','latex','FontSize',10); grid;
xlim([1 N]);
toc;

% Save the results
save CT_XPA.mat