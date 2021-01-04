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
maxitK = 300;     % maximum number of iterations in the outer loop
crit = 1e-6;      % criterion for the inner loop
critK = 1e-5;     % criterion for the outer loop
Delta = 1000;     % delta in HJB algorithm
% damp = 0.001;     % relaxation parameter for the deterministic steady state
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

% Set aggregate productivity grid
global Kmax Kmin intK intKK gridK dK
global Zmax Zmin Zmean intZ zmu zsigma gridZ dZ ddZ

% COMMENT : What is the optimum number of grids?
intK = 3;
intZ = 3; 

%Zmax = 2*sigma; Zmin = -2*sigma; intZ = 3; Zmean = 0;
Zmax = 2*sigma; Zmin = -2*sigma; Zmean = 0;
gridZ = linspace(Zmin,Zmax,intZ)'; dZ = (Zmax - Zmin)/(intZ - 1); ddZ = dZ^2;
gridZ((intZ+1)/2,1) = Zmean;

% Aggregate Shock Process
zmu = mu.*(Zmean - gridZ); 
zsigma = sigma^2.*ones(intZ,1);

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
Zshocks = randn(Stime,1);
% unused
% vtime = linspace(0,T,N); % ?

% Criteria for the outer loop
% epsmin = 1e-6; 
% epsilon = 1e+6; 
% iteration = 1; 