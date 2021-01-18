% Figure 1
clear all;

epsprint = 1;

load ./XPA/CT_XPA_sigma0.0070_kub0.20_klb0.20_zm4.00_intK3_intZ3.mat;
XPA_LOM = Kdot;

load ./KS/CT_KS_sigma0.0070_kub0.20_klb0.20_zm4.00_intK3_intZ3.mat;
KS_LOM = Kdot;

figure;
for i = 1:3
    subplot(2,3,i);
    plot(gridK,KS_LOM(:,i),'rx--','LineWidth',2);
    hold on;
    plot(gridK,XPA_LOM(:,i),'bo-','LineWidth',2);
    xlabel('$K$', 'interpreter','latex','FontSize',10); ylabel('$K''$', 'interpreter','latex','FontSize',10);
    title(sprintf('$Z=%1.4f$',gridZ(i)), 'interpreter','latex','FontSize',10);
    if (i==3); legend('KS','XPA','Location','SouthEast', 'interpreter','latex','FontSize',10); end;
    xlim([gridK(1) gridK(end)]); ylim([-0.4 0.4]);
%    xticks([34 35 36 37]);
end

if (epsprint); print -depsc2 Figure1.eps; end;

% Figure 2
load ./REITER/CT_REITER_sigma0.0070.mat;
REITER_K = KKpath';

% Plotting a figure comparing the path of aggregate capital for XPA, KS and REITER
figure;
plot(XPA_K(2:end,1),'b-','LineWidth',1);
hold on
plot(KS_K(2:end,1),'r:','LineWidth',1);
hold on
plot(REITER_K(1:end-2)','k--','LineWidth',1); 
title('Simulaton Path : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('XPA', 'KS', 'REITER', 'Location','northwest','interpreter','latex'); grid;

if (epsprint); print -depsc2 Figure2-1.eps; end;

% Plotting a figure comparing the path of aggregate capital for XPA, KS and REITER within a specific period
figure; 
plot(XPA_K(2:end,1),'b-','LineWidth',1);
hold on
plot(KS_K(2:end,1),'r:','LineWidth',1);
hold on
plot(REITER_K(1:end-2)','k--','LineWidth',1); 
title('Simulaton Path : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('XPA', 'KS', 'REITER', 'Location','northeast','interpreter','latex');
xlim([400 600]); grid;

if (epsprint); print -depsc2 Figure2-2.eps; end;