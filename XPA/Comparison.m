%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the comparison of Den Haan Error in the case of high TFP Volatility
% Masakazu Emoto @ Kobe univerisity 2020/07/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('XPA_DH.mat')
load('KS_DH.mat')
load('REITER_DH.mat')

figure(1)
plot(XPA_DH(:,1),'bo-','LineWidth',1);
hold on
plot(KS_DH(:,1),'ro-','LineWidth',1);
hold on
plot(REITER_DH(:,1),'ko-','LineWidth',1);
title('Den Haan Error', 'interpreter','latex','FontSize',10);
xlabel('Volatility of TFP : $\sigma$', 'interpreter','latex','FontSize',10);
ylabel('Den Haan Error : $\varepsilon^{MAX}$', 'interpreter','latex','FontSize',10);
legend('XPA', 'KS', 'REITER', 'Location','northwest','interpreter','latex');
xticks([1 2 3 4 5 6 7 8 9]);
xticklabels({'1','1.5','2','2.5','3','3.5','4','4.5','5'}); grid;