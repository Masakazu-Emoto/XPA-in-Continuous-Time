% Figure 1
clear all;
load ./XPA/CT_XPA_sigma0.0070_kub0.20_klb0.20_intK3_intZ3.mat;
XPA_LOM = Kdot;

load ./KS/CT_KS_sigma0.0070_kub0.20_klb0.20_intK3_intZ3.mat;
KS_LOM = Kdot;

figure;
for i = 1:3
    subplot(1,3,i);
    plot(gridK,KS_LOM(:,i),'rx--','LineWidth',2);
    hold on;
    plot(gridK,XPA_LOM(:,i),'bo-','LineWidth',2);
    xlabel('K','FontWeight','Normal'); ylabel('K''','FontWeight','Normal');
    title(sprintf('Z=%1.4f',gridZ(i)),'FontWeight','Normal');
    if (i==3); legend('KS','XPA','Location','SouthEast','FontWeight','Normal'); end;
%     xlim([gridK(1) gridK(end)]); ylim([-0.15 0.15]);
    xticks([34 35 36 37]);
end

%print -depsc2 LOM_210111.eps

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