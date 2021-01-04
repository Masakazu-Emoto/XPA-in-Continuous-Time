% main_KS.m : This code solves the Krusell and Smith model in continuous time by KS Algorithm

% Masakazu Emoto and Takeki Sunakawa (2021)
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time"
%
% Reference : Fernandez-Villaverde, J., S. Hurtado and G. Nuno (2019, FVHN hereafter)
% "Solving the Krusell-Smith (1998) model"
%
% The original code is downloaded from https://github.com/gregkaplan/phact

% Author : Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Revised by Takeki Sunakawa 2021/01/05
% E-mail address : masakazu.emoto@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALGORITHM
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
% Step 5 : Plot relevant graphs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; format long;
clc; tic;

%% -------------------------------------------------- %
%% Step 0 : Set Parameters
% -------------------------------------------------- %
global gamma rho alpha delta la intx x mu sigma com tau LAve 
global maxit maxitK crit critK Delta damp relax

% Preference
gamma = 1;        % Coefficient of relative risk aversion
rho = 0.01;       % Discount rate

% Production
alpha = 0.36;     % Capital share
delta = 0.025;    % Depriciate rate

% Idiosyncratic shock for labor productivity
intx = 2;
x1 = 0; x2 = 1;
% Level of labor producticvity 
% x(1): unemployed, x(2) : employed
x = [x1, x2]; 

% Transition probability for labor productivity 
% lambda1 : unemployment -> employment
% lambda2 : employment -> unemployment
lambda1 = 0.5;                                         
lambda2 = (lambda1 / (x(2) * 0.93 - x(1))*(x(2) - x(2) * 0.93));
la = [lambda1, lambda2];

% Aggregate shock for TFP (Ornstein-Uhlenbeck process)
mu = 0.25;        % Mean
sigma = 0.05;     % Variance 

% Tax system (Ahn et al.(2018))
com = 0.15; 
tau = (com / x(2)) * (la(2) / la(1));

% Average labour supply
LAve = (la(1) * x(2) + la(2) * x(1)) / (la(1) + la(2));

% PARAMETERS : We refer FVHN for the convergence criterion
maxit  = 100;     % maximum number of iterations in the inner loop
maxitK = 100;     % maximum number of iterations in the outer loop
crit = 1e-6;      % criterion for the inner loop
critK = 1e-5;     % criterion for the outer loop
Delta = 1000;     % delta in HJB algorithm
damp = 0.001;     % relaxation parameter for the deterministic steady state
relax = 0.7;      % relaxation parameter for the law of motion
relax_dot = 0.3;  % Initial weight in the relaxation algorithm for PLM convergence
relax1 = 0.9;     % reduction of the weight in the relaxation algorithm : wePLM = wePLM*wePLM1+wePLM2
relax2 = 0.005;   

% -------------------------------------------------- %
% Set grid
% -------------------------------------------------- %
% Individual wealth grid
global inta amin amax grida da aa aaa xx xxx Aswitch

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

%% -------------------------------------------------- %
%% Step 1 : Solve for the deterministic steady state
% -------------------------------------------------- %

% Calculate deterministic steady state
disp('Calcutating deterinistic steady state')
[rds, wds, Kds, Ads, uds, cds, pds, ids, Vds, gds] = steadystate();
toc;

disp('Deterministic steady state')
disp(Kds) 
disp(sum(gds'*grida*da))

% Set aggregate wealth and productivity grid, based on the deterministic
% steady state
global Kmax Kmin intK intKK gridK dK
global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ

% COMMENT : What is the optimum number of grids?
intK = 3;
intZ = 3; 

% Set the grid around the deterministic steady state of K
if sigma >= 0.03
    Kmax = 1.15*Kds; Kmin = 0.85*Kds; %intK = 3;
else
    Kmax = 1.05*Kds; Kmin = 0.95*Kds; %intK = 3;
end
gridK = linspace(Kmin,Kmax,intK)'; dK = (Kmax - Kmin)/(intK - 1);

%Zmax = 2*sigma; Zmin = -2*sigma; intZ = 3; Zmean = 0;
Zmax = 2*sigma; Zmin = -2*sigma; Zmean = 0;
gridZ = linspace(Zmin,Zmax,intZ)'; dZ = (Zmax - Zmin)/(intZ - 1); ddZ = dZ^2;
gridZ((intZ+1)/2,1) = Zmean;

% Aggregate Shock Process
zmu = mu.*(Zmean - gridZ); 
zsigma = sigma^2.*ones(intZ,1);

%% -------------------------------------------------- %
%% Step 2 : KS algorithm
% -------------------------------------------------- %
% Initial guess for the forecasting rule: the same as in FVHN
% What are KKdot and PLM_visits???
Kdot = zeros(intK,intZ); 
%KKdot = zeros(intKK, intZ);
%PLM_visits = zeros(intKK, intZ).*NaN; %rng(100);  

% PLM = zeros(inta,intx,intK,intZ);
% for ia = 1:inta
%     for ix = 1:intx
%         PLM(ia,ix,:,:) = Kdot(:,:);
%     end
% end

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

% Initial guess of vss, r, and w for the inner loop (used in inner.m)
r = alpha*quadK.^(alpha - 1).*LAve^(1 - alpha).*exp(quadZ) - delta;
w = (1 - alpha)*quadK.^(alpha).*LAve^(-alpha).*exp(quadZ);
if gamma == 1
    vss = log((w.*(1 - tau).*quadx + w.*com.*(1 - quadx) + r.*quada))/rho;
else
    vss = (w.*(1 - tau).*quadx + w.*com.*(1 - quadx) + r.*quada).^(1-gamma)/(1-gamma)/rho;
end

% Shock for Aggregate productivity
global T N Stime Dtime vtime dT
% the below is used in simulation.m
T = 500;  % years
N = 2000; % periods for simulation
dT = T/N; % interval of time, = 0.25
% the below is used in fokker_planck.m
rng(100); 
Stime = 1000; % the length of simulation 
Dtime = 500; % the length of simulation discarded to estimate the PLM
muini = gds; % stationary distribution from the deterministic steady state
Zshocks = randn(Stime,1);
% unused
vtime = linspace(0,T,N); % ?

% Inner loop and outer loop
disp('Calculating the inner and outer loops by KS algorithm')
disp(' ')
% Criteria for the algorithm
epsmin = 1e-6; 
epsilon = 1e+6; 
iteration = 1; 

while (epsilon > epsmin)
    %% -------------------------------------------------- %
    %% Step 2-1 : Inner Loop, Calculate the policy function by taking the forecasting rule (perceived law of motion) as given 
    %% -------------------------------------------------- %
    [A1, Bss, A3, vss, cs, ps] = inner_org(Kdot, vss, r, w);
    disp('Finished solving the HJB Equation')
    
    %% -------------------------------------------------- %
    %% Step 2-2 : Outer Loop (1), Simulate the path of aggregate capital
    %% -------------------------------------------------- %
%    [Ksim, Zsim, Zdown, Zup, Kdown, Kup] = fokker_planck(Zshocks, muini, Bss);
    [Ksim, Zsim, Kdown, Kup, Zdown, Zup] = fokker_planck_v1(Zshocks, muini, Bss);
    disp('Finished simulating aggregate capital')
    
    %% -------------------------------------------------- %
    %% Step 2-3 : Outer Loop (2), Solve for the forecasting rule by linear
    % regression of simulated data
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
    B = (X'*X)\(X'*Y); % Coefficient for Regression
    Y_LR = X * B;
    Y_Stad = (mean(Y - Y_LR).^2).^0.5;
    Y_R2 = 1 - (sum((Y - Y_LR).^2))/(sum((Y - mean(Y)).^2));
    
    % TS: Why do we use quadK and quadZ for updating Kdot???
    % what are these???
    X1mm = squeeze(quadK(1,1,:,:)) ;
    X2mm = squeeze(quadZ(1,1,:,:)) ;

    X1m = reshape(X1mm,[intK*intZ,1]);
    X2m = reshape(X2mm,[intK*intZ,1]);
    X0m = ones(size(X1m));
    
    X1m=log(X1m); 
    %X3m = X1m.*X2m;
    X_LR = [X0m X1m X2m];
    %X_LR = [X0m X1m X2m X3m];
    
    Kdotnew = X_LR*B; % B is obtained by the regression
    Kdotnew = reshape(Kdotnew, size(Kdot)); % Same value ! ????????
    
    if isreal(Kdotnew) ~= 1
        disp('')
        disp('The matrix for aggregate dynamics is not real')
        disp('')
        break
    end
    
    epsilon = max(max(abs(Kdot - Kdotnew)));
    iteration = iteration + 1;
    
    if epsilon > epsmin
%        disp('Calculating the law of motion ...')
%        disp([Y_R2, epsilon])
        disp([iteration, epsilon, Y_R2])
    end
    
    % Update law of motion
    % TS: Why do we do this?
    Kdot = (1 - relax_dot) * Kdot + relax_dot * Kdotnew;
    relax_dot = relax_dot * relax1 + relax2;
    
    % TS: PLM is not used inside the loop
%     for ia = 1:inta
%         for ix = 1:intx
%             PLM(ia,ix,:,:) = Kdot(:,:);
%         end
%     end
end

toc;

% extract Kdot to PLM
PLM = zeros(inta,intx,intK,intZ);
for ia = 1:inta
    for ix = 1:intx
        PLM(ia,ix,:,:) = Kdot(:,:);
    end
end

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
N = 10000; 
muini = gds; 
Zsim = zeros(N,1); 
rng(100);  
shock = randn(N,1); 
shock(1,1) = 0; % why ???
mmu = -1 + mu;
for time = 1:N-1
    if time == 1
        Zsim(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT)); % Ahn et al
    else
        Zsim(time+1) = (1 - mmu * dT)^(-1) * (Zsim(time) + sigma * shock(time) * sqrt(dT)); % Ahn et al
    end
end
simTFP = exp(Zsim);
[sim_mu, KS_K, KS_KK] =  simulate(Zsim, muini, Bss, Kdotnew);

% KS_KK is simulated results using the dynamics of aggregate capital
% KS_K is simulated results using the dynamics of aggregate capital and the HJB equation
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
%% Step 5 : Plot relevant graphs
%% -------------------------------------------------- %

% Plotting law of motion
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

% Ploting Sparse Matrix for value function
figure(3)
spy(A3,'b');    

% Plot Transition Dynamics DSS to SSS
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