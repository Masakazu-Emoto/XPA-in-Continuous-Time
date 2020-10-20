%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time."
% By Masakazu Emoto and Takeki Sunakawa
% This code is intended to plot a figure comparing the respective capital paths of XPA, KS, and REITER from the simulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading Data
load('XPApath.mat')
load('KSpath.mat')
load('REITERpath.mat')

% Plotting a figure comparing the path of aggregate capital for XPA, KS and REITER
figure(1)
plot(XPApath(2:end,1),'b-','LineWidth',1);
hold on
plot(KSpath(2:end,1),'r:','LineWidth',1);
hold on
plot(REITERpath(:,1:end-2)','k--','LineWidth',1); 
title('Simulaton Path : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('XPA', 'KS', 'REITER', 'Location','northwest','interpreter','latex'); grid;

% Plotting a figure comparing the path of aggregate capital for XPA, KS and REITER within a specific period
figure(2) 
plot(XPApath(2:end,1),'b-','LineWidth',1);
hold on
plot(KSpath(2:end,1),'r:','LineWidth',1);
hold on
plot(REITERpath(:,1:end-2)','k--','LineWidth',1); 
title('Simulaton Path : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('XPA', 'KS', 'REITER', 'Location','northeast','interpreter','latex');
xlim([400 600]);  grid;
