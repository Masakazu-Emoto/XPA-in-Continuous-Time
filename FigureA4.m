% Figure A4
clear all;

epsprint = 1;

load ./XPA/CT_XPA_sigma0.0070_kub0.20_klb0.20_intK3_intZ3.mat;

figure
plot(XPA_K,'b-','LineWidth',1);
hold on
plot(XPA_KK,'r-','LineWidth',1);
hold on
title('Simulaton Path : XPA : $\sigma$ = 0.007', 'interpreter','latex','FontSize',10);
xlabel('Simulation time : $T$', 'interpreter','latex','FontSize',10);
ylabel('Capital : $K$', 'interpreter','latex','FontSize',10);
legend('$K^*_{t}$', '$\tilde{K}_{t}$','Location','northwest','interpreter','latex'); grid;
xlim([1 N]);

if (epsprint); print -depsc2 FigureA4.eps; end;