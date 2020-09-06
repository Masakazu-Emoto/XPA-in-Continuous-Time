
load('XPApath.mat')
load('KSpath.mat')
load('REITERpath.mat')

figure(1)
plot(XPApath(2:end,1),'b-','LineWidth',1);
hold on
plot(KSpath(2:end,1),'r-','LineWidth',1);
hold on
plot(REITERpath(:,1:end-2)','k-','LineWidth',1); 
title('Simulaton Path :: $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('XPA', 'KS', 'REITER', 'Location','northwest','interpreter','latex'); grid;

figure(2)
plot(XPApath(2:end,1),'b-','LineWidth',1);
hold on
plot(KSpath(2:end,1),'r-','LineWidth',1);
hold on
plot(REITERpath(:,1:end-2)','k-','LineWidth',1); 
title('Simulaton Path :: $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('XPA', 'KS', 'REITER', 'Location','northeast','interpreter','latex');
xlim([400 600]);  grid;
