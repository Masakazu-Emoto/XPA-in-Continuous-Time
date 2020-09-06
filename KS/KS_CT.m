%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Krusell-Smith Model (1998) in Continuous Time by KS Algorithm
% KS algorithm is refered by Villaverde et al.(2019 NBER)
% Masakazu EMOTO @ Kobe univerisity 2020/07/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm 
% Phase 1 : Calculate deterministic steady state
% Phase 2 : Calculate the bias correction terms
% Phase 3-1 : Inner Loop Taking as given the forecasting rule, calculating policy function
% Phase 3-2 : Outer Loop Calculation Law of Motion
% Phase 4 : Calculate stochastic steady state
% Phase 5 : Simulation and Calculate Den Haan Error
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
sigma = 0.007; % Variance 

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
relax_dot   = 0.3;  % Initial weigth in the relaxation algorithm for PLM convergence
relax1   = 0.9;       % reduction of the relaxation algorithm: wePLM = wePLM*wePLM1+wePLM2
relax2   = 0.005;    % reduction of the relaxation algorithm

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
[rds, wds, Kds, Ads, uds, cds, pds, ids, Vds, gds] = steadystate();
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
global Kmax Kmin intK intKK gridK dK
global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ

if sigma >= 0.03
    Kmax = 1.15*Kds; Kmin = 0.85*Kds; intK = 3;
else
    Kmax = 1.05*Kds; Kmin = 0.95*Kds; intK = 3;
end
gridK = linspace(Kmin,Kmax,intK)'; dK = (Kmax - Kmin)/(intK - 1);

Zmax = 2*sigma; Zmin = -2*sigma; intZ = 3; Zmean = 0;
gridZ = linspace(Zmin,Zmax,intZ)'; dZ = (Zmax - Zmin)/(intZ - 1); ddZ = dZ^2;
gridZ((intZ+1)/2,1) = Zmean;

% Aggregate Shock Process
zmu = mu.*(Zmean - gridZ); 
zsigma = sigma^2.*ones(intZ,1);

% -------------------------------------------------- %
% Phase 3 : KS algorithm
% -------------------------------------------------- %
% Initial guess for law of motion for Krusell and Smith (1998)
% Initial guess is the same as Villaverde et al.(2019)
Kdot = zeros(intK,intZ); KKdot = zeros(intKK, intZ);
PLM_visits = zeros(intKK, intZ).*NaN; %rng(100);  

for ia = 1:inta
    for ix = 1:intx
        PLM(ia,ix,:,:) = Kdot(:,:);
    end
end

% Inner loop and outer loop
disp('Calcutating inner loop and outer loop by KS algorithm')
% epsmin = 0.00050; epsilon = 1e+6; iteration = 1; 
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

% Initial guess for KS Algorithm
r = alpha*quadK.^(alpha - 1).*LAve^(1 - alpha).*exp(quadZ) - delta;
w = (1 - alpha)*quadK.^(alpha).*LAve^(-alpha).*exp(quadZ);
if gamma == 1
    vss = log((w.*(1 - tau).*quadx + w.*com.*(1 - quadx) + r.*quada))/rho;
else
    vss = (w.*(1 - tau).*quadx + w.*com.*(1 - quadx) + r.*quada).^(1-gamma)/(1-gamma)/rho;
end

% Shock for Aggregate productivity
global T N Stime Dtime vtime dT
T = 500; N = 2000; dT = T/N; rng(100); 
Stime =1000; Dtime = 500; vtime = linspace(0,T,N);
muini = gds; Zshocks = randn(Stime,1);

% Calculation for Forecasting Rules
while (epsilon > epsmin)
    PLM_visits_=PLM_visits;

    % -------------------------------------------------- %
    % Phase 3-1 : Calculate policy function : Inner Loop
    % -------------------------------------------------- %
    [Ass, WW, vss, cs, ps] = inner(Kdot, vss, iteration, r, w);
    disp('Finished calculation HJB Equation')
    
    % -------------------------------------------------- %
    % Phase 3-2 : Calculate law of motion : Outer Loop:Simulation
    % -------------------------------------------------- %
    [Ksim, Zsim, Zdown, Zup, Kdown, Kup]  = fokker_planck(Zshocks, muini, Ass);
    disp('Finished calculation Forecasting Rule')
    
    % -------------------------------------------------- %
    % Phase 3-3: Calculate law of motion : Linear Regression
    % -------------------------------------------------- %
    % Calculation linear regression
    Y = (Ksim(Dtime+1:Stime) - Ksim(Dtime:Stime-1))/dT; % dK(t): Growth rate of K at time t
    X0 = ones(Stime-Dtime,1);         % Constant
    X1 = Ksim(Dtime:Stime-1);         % K(t-1)
    X2 = Zsim(Dtime:Stime-1);    % Z(t-1)
    X1 = log(X1);
    X = [X0 X1 X2];
    
    B = (X'*X)^-1*X'*Y;   % Coefficient for Regression
    Y_LR = X * B;
    Y_Stad = (mean(Y - Y_LR).^2).^0.5;
    Y_R2 = 1 - (sum((Y - Y_LR).^2))/(sum((Y - mean(Y)).^2));
    
    X1mm = squeeze(quadK(1,1,:,:)) ;
    X2mm = squeeze(quadZ(1,1,:,:)) ;

    X1m = reshape(X1mm,[intK*intZ,1]);
    X2m = reshape(X2mm,[intK*intZ,1]);
    X0m = ones(size(X1m));
    
    X1m=log(X1m); 
    X_LR = [X0m X1m X2m];
    Kdotnew = X_LR*B;
    Kdotnew = reshape(Kdotnew, size(Kdot)); % Same value !
    
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
        disp([Y_R2, epsilon])
    end
    
    % Update law of motion
    Kdot = (1 - relax_dot) * Kdot + relax_dot * Kdotnew;
    relax_dot = relax_dot * relax1 + relax2;
    
    for ia = 1:inta
        for ix = 1:intx
            PLM(ia,ix,:,:) = Kdot(:,:);
        end
    end
end

toc;

% -------------------------------------------------- %
% Phase 4 : Calculate stochastic steady state and Den Haan Error
% -------------------------------------------------- %
disp('Calculating stochastic steady state')
muini = gds; simZ = ones(N,1).*Zmean;
simTFP = exp(simZ);
[simvalue, sim_mu, Kss, KKss] = stochastic_steady(simZ, muini, Ass, Kdotnew);
toc;

% -------------------------------------------------- %
% Phase 5 : Simulation and Calculate Den Haan Error
% -------------------------------------------------- %
disp('Simulating the model and Calculating Den Haan Error')
N = 10000; muini = gds; Zsim = zeros(N,1); rng(100);  
shock = randn(N,1); shock(1,1) = 0;
mmu = -1 + mu;
for time = 1:N-1
    if time == 1
        Zsim(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT)); % Ahn et al method
    else
        Zsim(time+1) = (1 - mmu * dT)^(-1) * (Zsim(time) + sigma * shock(time) * sqrt(dT)); % Ahn et al method
    end
end
irfTFP = exp(Zsim);
[testvalue, test_mu, KS_K, KS_KK] =  simulate(Zsim, muini, Ass, Kdotnew);

% irf_KK : irf_KK is simulated results using the dynamics of aggregate capital
% irf_K : irf_K is simulated results using the dynamics of aggregate capital and the HJB equation
Drop = 1000; 
DH_Error = 100.0 * max(abs(log(KS_KK(Drop+1:end)) - log(KS_K(Drop+1:end))));
DH_Mean = 100.0 * sum(abs(log(KS_KK(Drop+1:end)) - log(KS_K(Drop+1:end))))/(N - Drop);

disp('MAX Den Haan Error')
disp(DH_Error)
disp('MEAN Den Haan Error')
disp(DH_Mean)
toc;

disp('Stochastic steady state by KS Algorithm')
disp(Kss(end)) 

% Plot Forecasting rule
figure(1)
surf(gridZ,gridK,Kdot);
title('The perceived law of motion, PLM : KS', 'interpreter','latex','FontSize',14);
xlabel('shock ($Z$)', 'interpreter','latex','FontSize',14);
ylabel('capital ($K$)', 'interpreter','latex','FontSize',14);
xlim([Zmin Zmax]); ylim([Kmin Kmax]);

intKK = 16; gridKK = linspace(Kmin,Kmax,intKK)';
intZZ = 41; gridZZ = linspace(Zmin, Zmax,intZZ)';
LOM = interp1(gridK,Kdot,gridKK,'spline');
LOM = interp1(gridZ,LOM',gridZZ,'spline');

figure(2)
surf(gridZZ,gridKK,LOM');
title('Forecasting rule : KS : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Shock : $Z$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
zlabel('$\dot{K}$', 'interpreter','latex','FontSize',10);
xlim([Zmin Zmax]); ylim([Kmin Kmax]);

% Plot Sparse Matrix
figure(3)
spy(WW,'b');    

% Plot Transition Dynamics DSS to SSS
figure(4)
plot(Kss,'b-','LineWidth',1); grid;
hold on
plot(Kss(1,1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1); grid;
hold on
plot(max(size(Kss)),Kss(end,1),'bo','MarkerEdgeColor','b','MarkerFaceColor','b');
title('Transition Dynamics the DSS to the SSS : $\sigma$ = 0.07', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10); grid;

% Plot Simulaion Path
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
save CT_KS.mat