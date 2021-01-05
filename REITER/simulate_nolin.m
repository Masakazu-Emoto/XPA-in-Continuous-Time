function [linvalues, vtime, nonvalues] = simulate_nolin(g1,impact,T,N,shocks,method,blowup,subset)
%% simulate_nolin.m : This code modifies simulate.m so as to calculate the Den Haan Error
%    Line 72-82 : We calculate the case of the non-linear dynamics about the wealth distribution dynamics by Implict method
%    Line 100-110 : We calculate the case of the non-linear dynamics about the wealth distribution dynamics by Explict method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Given linear dynamics
%        dx = g1*x*dt + impact * dZ
%    simulates the values of x for T time period with N steps with
%    realization of dZ given by shocks
%
% by SeHyoun Ahn, March 2017
%
% REFERENCE: Ahn, SeHyoun, Greg Kaplan, Benjamin Moll, Thomas Winberry, and
%    Christian Wolf. "When Inequality Matters for Macro and Macro Matters
%    for Inequality."
%
% PARAMETERS:
%    dx = g1*x*dt + impact*dZ
%    g1 = (m x m) matrix for dynamics, usually an output of
%         <schur_solver.m>
%    impact = (m x n) matrix of impact of shocks on values, usually an
%             output of <schur_solver.m>
%    T = length of time simulation
%    N = number of steps taken for the simulation
%    shocks = (n x N) matrix of shock realization
%    method = {'implicit','explicit'} updating method for linear system
%    blowup = (optional) translation back to full system if the dynamics is
%             given for a reduced system
%    subset = (optional) subset out variables if only part of the system is
%             needed
%
% OUTPUTS:
%    values = (m x N) matrix of simulated values (subset of the values if
%             optional subset parameter is used)
%    vtime = (N x 1) vector of time values
%
% EXAMPLES:
%    A = diag(-linspace(1,11,10));     % simple stable dynamics
%    impact = linspace(1,11,10)';      % simple impact of shocks
%    T = 5;
%    N = 50;
%    shocks = zeros(1,N); shocks(1) = 1;
%    [vtime,values] = simulate(A,impact,T,N,shocks,'implicit');
%    plot(vtime,values);
%
%    For example with blowup and subset, see Krusell-Smith case example
%       available at < > 
%
% SYNTAX:
% [values,vtime] = simulate(g1,impact,T,N,shocks,method,blowup,subset)

global ggamma rrho ddelta aalpha ssigmaTFP rrhoTFP z lla mmu ttau I amin amax a da aa ...
	zz Aswitch rmin rmax r0 maxit crit Delta Ir crit_S IfSS IbSS I0SS aaa zzz varsSS zAvg nVars nEErrors ggSS
vtime = linspace(0,T,N);
dt = vtime(2)-vtime(1);

% Preallocation
[nvars,~] = size(g1);
linvalues = zeros(nvars,N); nonvalues = zeros(4*I,N+1);
if (method == 'implicit')     
    % if N is small, might be faster just to do backslash instead. To use
    % backslash method, just uncomment line 51/54 and comment line 50/53
    
    % Linearized Solution
    gg1 = inv(speye(size(g1)) - g1*dt);
    for n = 1 : N
        linvalues(:,n+1) = gg1*(linvalues(:,n) + (dt^(1/2))*impact*shocks(:,n));
    end
    
    % For Den Haan Error : Non-Linearized Solution
    vdot(:,1) = linvalues(1:2*I,1); gdot(:,1) = linvalues(2*I+1:4*I,1); TFP(1,1) = linvalues(4*I,1);
    nonvalues(:,1) = [vdot(:,1); gdot(1:end-1,1); TFP(1,1)];
    for n = 1 : N
        [A,g,g_end] = calcA(gdot,vdot,TFP,n);
        gnext = (speye(I*2) - A' * dt)\[g; g_end]; % Implicit Case
        gdot(:,n+1) = gnext - varsSS(2*I+1:4*I);
        vdot(:,n+1) = gg1(1:2*I,:)*(nonvalues(:,n) + (dt^(1/2))*impact*shocks(:,n));
        TFP(:,n+1) = gg1(4*I,:)*(nonvalues(:,n) + (dt^(1/2))*impact*shocks(:,n));
        nonvalues(:,n+1) = [vdot(:,n+1); gdot(1:end-1,n+1); TFP(1,n+1)];
    end
    linvalues = linvalues(:,2:end);
    nonvalues = nonvalues(:,2:end);
    if nargin > 6
        linvalues = blowup*linvalues;
        if nargin > 7
            linvalues = linvalues(subset,:);
        end
    end
    
elseif (method == 'explicit')

    % Linearized Solution    
    gg1 = speye(size(g1))+g1*dt;
    for n = 1 : N
        linvalues(:,n+1) = gg1*(linvalues(:,n) + (dt^(1/2))*impact*shocks(:,n));	
    end
    
    % For Den Haan Error : Non-Linearized Solution
    % NOTE: 1 instead of 2?
    vdot(:,1) = linvalues(1:2*I,2); gdot(:,1) = linvalues(2*I+1:4*I,2); TFP(1,1) = linvalues(4*I,2);
    nonvalues(:,1) = [vdot(:,1); gdot(1:end-1,1); TFP(1,1)];
    for n = 1 : N
        [A,g,g_end] = calcA(gdot,vdot,TFP,n);
        gnext = [g; g_end] + A' * [g; g_end]*dt; % Explicit Case
        gdot(:,n+1) = gnext - varsSS(2*I+1:4*I);
        vdot(:,n+1) = gg1(1:2*I,:)*(nonvalues(:,n) + (dt^(1/2))*impact*shocks(:,n));
        TFP(:,n+1) = gg1(4*I,:)*(nonvalues(:,n) + (dt^(1/2))*impact*shocks(:,n));
        nonvalues(:,n+1) = [vdot(:,n+1); gdot(1:end-1,n+1); TFP(1,n+1)];
    end
    linvalues = linvalues(:,2:end);
    nonvalues = nonvalues(:,2:end);
    if nargin > 6
        linvalues = blowup*linvalues;
        if nargin > 7
            linvalues = linvalues(subset,:);
        end
    end
else
    error('<simulate>: Unknown method');
end

end

function [A,g,g_end] = calcA(gdot,vdot,TFP,n)
% calculate the transition matrix

    global ggamma rrho ddelta aalpha ssigmaTFP rrhoTFP z lla mmu ttau I amin amax a da aa ...
    zz Aswitch rmin rmax r0 maxit crit Delta Ir crit_S IfSS IbSS I0SS aaa zzz varsSS zAvg nVars nEErrors ggSS
    
    g = gdot(1:end-1,n) + varsSS(2*I+1:4*I-1);	% vector of distribution, removing last point		
    g_end = 1/da-sum(g);

    V = vdot(1:2*I,n) + varsSS(1:2*I);
    V = reshape(V,I,2);

    T = TFP(1,n);
    K = sum(aaa .* [g;g_end] * da); 
    r = exp(T) * aalpha * (K ^ (aalpha - 1)) * (zAvg ^ (1 - aalpha)) - ddelta;
    w = exp(T) * (1 - aalpha) * (K ^ aalpha) * (zAvg ^ (-aalpha)); 

    % Compute forward difference
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = (w*((1 - ttau) * z + mmu * (1 - z)) + r.*amax).^(-ggamma); %will never be used, but impose state constraint a<=amax just in case

    % Compute backward difference
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = (w*((1 - ttau) * z + mmu * (1 - z)) + r.*amin).^(-ggamma); %state constraint boundary condition

    % Compute consumption and savings with forward difference
    cf = dVf.^(-1/ggamma);
    ssf = w*((1 - ttau) * zz + mmu * (1 - zz)) + r.*aa - cf;

    % Compute consumption and savings with backward difference
    cb = dVb.^(-1/ggamma);
    ssb = w*((1 - ttau) * zz + mmu * (1 - zz)) + r.*aa - cb;

    % Compute consumption and derivative of value function for no drift
    c0 = w*((1 - ttau) * zz + mmu * (1 - zz)) + r.*aa;
    dV0 = c0.^(-ggamma);

    % Compute upwind differences    
    If = ssf > 0;       % positive drift --> forward difference
    Ib = ssb < 0;       % negative drift --> backward difference
    I0 = (1-If-Ib);     % no drift
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0;
    c = dV_Upwind.^(-1/ggamma);

    % Construct matrix for updating implicit scheme
    X = -min(ssb,0)/da;
    Y = -max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;

    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    
end
