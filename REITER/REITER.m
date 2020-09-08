%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Krusell-Smith Model (1998) in Continuous Time by Ahn-Reiter Algorithm
% This code modifies Ahn's GitHub code enables to calculate the Den Haan Error
% Masakazu EMOTO @ Kobe univerisity 2020/07/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all; close all;
tstart = tic;

%% Set options for this example run
% initialize shocks for simulation
T = 2500; N = 10000; rng(100)
vAggregateShock = randn(1,N+1);    

%% Step 0: Set Parameters
% The script sets up parameters relevant for the model
set_parameters;

%% Step 1: Solve for Steady State
% Non-stochastic steady state can be found using any methods. In
%    particular, example codes can be found at
%    <<http://www.princeton.edu/%7Emoll/HACTproject.htm>>.

tStart = tic;
fprintf('Computing steady state...\n')
global IfSS IbSS I0SS varsSS A ggSS

[rSS,wSS,KSS,ASS,uSS,cSS,VSS,gSS,dVUSS,dVfSS,dVbSS,IfSS,IbSS,I0SS] = compute_steady_state();

fprintf('Time to compute steady state: %.3g seconds\n\n\n',toc(tStart));

% Store steady state values
varsSS = zeros(nVars,1);
varsSS(1:2*I,1) = reshape(VSS,2*I,1);
ggSS = reshape(gSS,2*I,1);
varsSS(2*I+1:4*I-1,1) = ggSS(1:2*I-1);
varsSS(4*I,1) = 0;
varsSS(4*I+1,1) = KSS;
varsSS(4*I+2,1) = rSS;
varsSS(4*I+3,1) = wSS;
varsSS(4*I+4,1) = (KSS ^ aalpha) * (zAvg ^ (1 - aalpha));
CSS = sum(cSS(:) .* gSS(:) * da);
varsSS(4*I+5,1) = CSS;
varsSS(4*I+6,1) = ddelta * KSS;

%% Step 2: Linearize Model Equations
% For computing derivatives, the codes written for solving for the
%    steady-state can be used almost verbatim using automatic
%    differentiation toolbox as long as only the functions supported by
%    automatic differentation are used. For list of supported functions and
%    documentation of relevant syntax check <<https://github.com/sehyoun/MATLABAutoDiff>>
fprintf('Taking derivatives of equilibrium conditions...\n')
t0 = tic;

% Prepare automatic differentiation
vars = zeros(2*nVars+nEErrors+1,1);
vars = myAD(vars);

% Evaluate derivatives
[derivativesIntermediate,VEError] = equilibrium_conditions(vars);

% Extract out derivative values
derivs = getderivs(derivativesIntermediate);
Vderivs = getderivs(VEError);
tDerivs = toc(t0);
fprintf('...Done!\n')
fprintf('Time to compute derivatives: %2.4f seconds\n\n\n',tDerivs)
if tDerivs > 1
    warning('If you compile mex files for automatics differentiation, matrix vector multiplication will be slow');
    disp('Press any key to continue...');
    pause();
end

% Unpackage derivatives
mVarsDerivs = derivs(:,1:nVars);
mVarsDotDerivs = derivs(:,nVars+1:2*nVars);
mEErrorsDerivs = derivs(:,2*nVars+1:2*nVars+nEErrors);
mShocksDerivs = derivs(:,2*nVars+nEErrors+1);

%% Step 3: Solve Out Static Constraints

% rename derivatives to match notation in paper
g0 = mVarsDotDerivs;
g1 = -mVarsDerivs;
c = sparse(nVars,1);
psi = -mShocksDerivs;
pi = -mEErrorsDerivs;

t0 = tic;
fprintf('Model Reduction ...\n')

% Clean G0
[state_red,inv_state_red,g0,g1,c,pi,psi] = clean_G0_sparse(g0,g1,c,pi,psi);
n_g_red = n_g;
from_spline = speye(n_g_red + n_v);
to_spline = speye(n_g_red + n_v);
n_splined = n_v;

%% Step 4: Solve Linear System
t0 = tic;
fprintf('Solving reduced linear system...\n')

[G1,~,impact,eu,F] = schur_solver(g0,g1,c,psi,pi,1,1,1);

fprintf('...Done!\n')
fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
fprintf('Time to solve linear system: %2.4f seconds\n\n\n',toc(t0))
toc(tstart)
%% Step 5: Simulate the model and Calcuate the Den Haan Error
fprintf('Simulating Model...\n')
t0 = tic;

trans_mat = inv_state_red*from_spline;
[simulated,vTime,G,vvalues] = simulate(G1,impact,T,N,vAggregateShock,'implicit',inv_state_red);

Gwealth = G(1:2*I-1,:) + ggSS(1:end-1);
G_End = 1/da - sum(Gwealth);
Gpath = sum(aaa .* [Gwealth; G_End] * da);

big_sim = simulated + varsSS;
big_value = vvalues(:,2:end) + varsSS(1:4*I,1);
big_TFP = vvalues(end,2:end);
fprintf('...Done!\n')
fprintf('Time to simulate model: %2.4f seconds\n\n\n',toc(t0))

% Add state-states back in to get values in levels
varsSS_small = varsSS(4*I:4*I+6,1);
vAggregateTFP = simulated(400,:) + varsSS_small(1);
Kpath = simulated(401,:) + varsSS_small(2);
vAggregateOutput = simulated(404,:) + varsSS_small(5);
vAggregateConsumption = simulated(405,:) + varsSS_small(6);
vAggregateInvestment = simulated(406,:) + varsSS_small(7);

% Compute log differences for plotting
vAggregateTFP_reduced = vAggregateTFP;
vAggregateOutput_reduced = log(vAggregateOutput) - log(varsSS_small(5));
vAggregateConsumption_reduced = log(vAggregateConsumption) - log(varsSS_small(6));
vAggregateInvestment_reduced = log(vAggregateInvestment) - log(varsSS_small(7));

%% (optional) Step 7: Plot relevant values
Drop = 1000;
DH_Error = 100.0*max(abs(log(Kpath(1001:end-1)) - log(Gpath(1002:end))))
DH_Mean = 100.0 * sum(abs(log(Kpath(1001:end-1)) - log(Gpath(1002:end))))/(N - Drop)

% Simulation path Linear and Non-linear
figure(2)
plot(Gpath(2:end),'b-','LineWidth',1.5);
hold on
plot(Kpath(1:end-1),'r-','LineWidth',1.5);
title('Simulaton Path : REITER : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10); grid;
legend('Non-Linear','Linear','Location','northwest','interpreter','latex');

toc(tstart)