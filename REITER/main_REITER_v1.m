%% main_REITER_v1.m : These are programs for solving Krusell-Smith model in continuous time by the REITER algorithm as described in
%
% Masakazu Emoto and Takeki Sunakawa (2021)
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time"
%
% Reference : Ahn et al. (2018, NBER Macroeconomics Annual 32.1, 1-75)
% "When inequality matters for macro and macro matters for inequality"
%
% The original code is downloaded from
% https://github.com/gregkaplan/phact
%
% Author : Masakazu EMOTO @ Kobe univerisity 2020/10/22
% Revised by Takeki Sunakawa 2021/01/05
% E-mail address : masakazu.emoto@gmail.com
%
% Uses : autodiff toolbox, phact toolbox (including our simulate_nolin.m for the Den Haan Errors),
% set_parameters.m, compute_steady_state.m, equilibrium_conditions.m
% NOTE: We modify the original simulate.m in the phact toolbox to calculatethe Den Haan Errors.
%
%% Summary of the algorithm
%    Step 0 : Set Parameters
%    Step 1 : Solve for Deterministic Steady State
%    Step 2 : Linearize Model Equations
%    Step 3 : Solve Out Static Constraints
%    Step 4 : Solve Linear System
%    Step 5 : Simulate the model and Calcuate the Den Haan Error
%    Step 6 : Plot relevant values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all;
% close all;
tstart = tic;

diagnose = 0;
loadtemp = 1;
savetime = 0;

%% Set options for this example run
% initialize shocks for simulation
%% NOTE: This is used for calculating the Den Haan Errors with simulate_nolin.m
T = 2500; N = 10000; rng(100)
vAggregateShock = randn(1,N+1);

%% Step 0: Set Parameters
% The script sets up parameters relevant for the model
set_parameters;
if (loadtemp)
    load temp.mat ssigmaTFP;
    fprintf('std for Zshock = %1.4f\n',ssigmaTFP);
end

%% Step 1: Solve for Steady State
% Non-stochastic steady state can be found using any methods. In
%    particular, example codes can be found at
%    <<http://www.princeton.edu/%7Emoll/HACTproject.htm>>.

tStart = tic;
fprintf('Computing steady state...\n')
global IfSS IbSS I0SS varsSS A ggSS

[rSS,wSS,KSS,ASS,uSS,cSS,VSS,gSS,dVUSS,dVfSS,dVbSS,IfSS,IbSS,I0SS] = compute_steady_state();
if (diagnose); fprintf('Time to compute steady state: %.3g seconds\n\n\n',toc(tStart)); end;

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
t1 = tic;

% Prepare automatic differentiation
vars = zeros(2*nVars+nEErrors+1,1);
vars = myAD(vars);

% Evaluate derivatives
[derivativesIntermediate,VEError] = equilibrium_conditions(vars);

% Extract out derivative values
derivs = getderivs(derivativesIntermediate);
% Vderivs = getderivs(VEError); %% for what???
tDerivs = toc(t0);
fprintf('...Done!\n')
if (diagnose); fprintf('Time to compute derivatives: %2.4f seconds\n\n\n',tDerivs); end;
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

% t0 = tic;
% fprintf('Model Reduction ...\n')

% Clean G0
%% NOTE: We set reduceDistribution = 0 and reduceV = 0 in the original code
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
if (diagnose)
    fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
    fprintf('Time to solve linear system: %2.4f seconds\n\n\n',toc(t0));
end
%if (diagnose); toc(tstart); end;
etime1 = toc(t1);
%fprintf('...Done!\n');
if (diagnose); fprintf('Time to solve model = %2.4f\n',etime1); end;

if (~savetime)
%% Step 5: Simulate the model and Calcuate the Den Haan Error
    
    fprintf('Simulating Model...\n')
    t0 = tic;

    trans_mat = inv_state_red*from_spline;
    %[simulated,vTime] = simulate(G1,impact,T,N,vAggregateShock,'implicit',trans_mat,4*I:4*I+6);
    %% ------
    %% NOTE: We also calculate "non-linear dynamics" to obtain the Den Haan Errors
    [lin_simulated,vTime,non_lin_simulated] = simulate_nolin(G1,impact,T,N,vAggregateShock,'implicit',inv_state_red);
    %[~,~,non_lin_simulated] = simulate_nolin(G1,impact,T,N,vAggregateShock,'implicit',inv_state_red);

    % Non-linear dynamics
    gg = non_lin_simulated(2*I+1:4*I-1,:) + ggSS(1:end-1);
    gg_End = 1/da - sum(gg);
    KKpath = sum(aaa .* [gg; gg_End] * da);
    %% ------

    fprintf('...Done!\n')
    if (diagnose); fprintf('Time to simulate model: %2.4f seconds\n\n\n',toc(t0)); end;

    % Add state-states back in to get values in levels
    varsSS_small = varsSS(4*I:4*I+6,1);
    %% NOTE: indices are changed
    % vAggregateTFP = lin_simulated(1,:) + varsSS_small(1);
    % vAggregateOutput = lin_simulated(5,:) + varsSS_small(5);
    % vAggregateConsumption = lin_simulated(6,:) + varsSS_small(6);
    % vAggregateInvestment = lin_simulated(7,:) + varsSS_small(7);
    %
    % Kpath = lin_simulated(2,:) + varsSS_small(2);

    vAggregateTFP = lin_simulated(400,:) + varsSS_small(1);
    vAggregateOutput = lin_simulated(404,:) + varsSS_small(5);
    vAggregateConsumption = lin_simulated(405,:) + varsSS_small(6);
    vAggregateInvestment = lin_simulated(406,:) + varsSS_small(7);

    Kpath = lin_simulated(401,:) + varsSS_small(2);

    % Compute log differences for plotting
    vAggregateTFP_reduced = vAggregateTFP;
    vAggregateOutput_reduced = log(vAggregateOutput) - log(varsSS_small(5));
    vAggregateConsumption_reduced = log(vAggregateConsumption) - log(varsSS_small(6));
    vAggregateInvestment_reduced = log(vAggregateInvestment) - log(varsSS_small(7));

    %% Step 6: Plot relevant values

    % KKpath is the sequence of simulated results using the non-linear dynamics
    % Kpath is the sequence of simulated results using the linearized dynamics
    Drop = 1000;
    DH_Max = 100.0 * max(abs(log(Kpath(1001:end-1)) - log(KKpath(1002:end))));
    DH_Mean = 100.0 * sum(abs(log(Kpath(1001:end-1)) - log(KKpath(1002:end))))/(N - Drop);

    % % Simulation path Linear and Non-linear
    % figure
    % plot(KKpath(2:end),'b--','LineWidth',1.5);
    % hold on
    % plot(Kpath(1:end-1),'r-','LineWidth',1.5);
    % title(['Simulaton Path : REITER : $\sigma =$', num2str(ssigmaTFP,'%.3f')], 'interpreter','latex','FontSize',10);
    % xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
    % ylabel('Capital : $K$', 'interpreter','latex','FontSize',10); grid;
    % legend('Non-Linear','Linear','Location','northwest','interpreter','latex');

    if (loadtemp)
        disp('done');
        disp(' ');
        if (1-rrhoTFP==0.25)
            eval(sprintf('save CT_REITER_sigma%1.4f.mat',ssigmaTFP));
        else % robustness for mu
            eval(sprintf('save CT_REITER_mu%1.2f_sigma%1.4f.mat',1-rrhoTFP,ssigmaTFP));
        end    
    end
    
end

if (savetime); save etime.mat etime1; end;

%toc(tstart)
