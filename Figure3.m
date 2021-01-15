% Figure 3
clear all;

epsprint = 0;

kub = 0.2;
klb = 0.2;
zm = 4;
intK = 3;
intZ = 3;
% vsig = [0.01:0.001:0.03]';
vsig = [0.01:0.005:0.05]';
%vsig = [0.01:0.01:0.05]';

DH_KS = zeros(size(vsig,1),2);
DH_XPA = zeros(size(vsig,1),2);
DH_REITER = zeros(size(vsig,1),2);
for i = 1:size(vsig,1)
    
    sigma = vsig(i);
    eval(sprintf('load ./KS/CT_KS_sigma%1.4f_kub%1.2f_klb%1.2f_zm%1.2f_intK%d_intZ%d.mat',sigma,kub,klb,zm,intK,intZ));
    DH_KS(i,1) = DH_Max;
    DH_KS(i,2) = DH_Mean;
    eval(sprintf('load ./XPA/CT_XPA_sigma%1.4f_kub%1.2f_klb%1.2f_zm%1.2f_intK%d_intZ%d.mat',sigma,kub,klb,2.5,intK,intZ));
    DH_XPA(i,1) = DH_Max;
    DH_XPA(i,2) = DH_Mean;
    eval(sprintf('load ./REITER/CT_REITER_sigma%1.4f.mat DH_Max DH_Mean',sigma));
    DH_REITER(i,1) = DH_Max;
    DH_REITER(i,2) = DH_Mean;
    
end

%  Plotting a figure comparing the values of Den Haan Error for XPA, KS and REITER
figure
plot(100*vsig,DH_XPA(:,1),'bo-','LineWidth',1);
hold on
plot(100*vsig,DH_KS(:,1),'ro:','LineWidth',1);
plot(100*vsig,DH_REITER(:,1),'ko--','LineWidth',1);
title('Den Haan Error', 'interpreter','latex','FontSize',10);
xlabel('Volatility of TFP : $\sigma$ (\%)', 'interpreter','latex','FontSize',10);
ylabel('Den Haan Error : $\varepsilon^{MAX}$ (\%)', 'interpreter','latex','FontSize',10);
legend('XPA', 'KS', 'REITER', 'Location','northwest','interpreter','latex');
% xticks([1 2 3 4 5 6 7 8 9]);
% xticklabels({'1','1.5','2','2.5','3','3.5','4','4.5','5'}); grid;

if (epsprint); print -depsc2 Figure3.eps; end;
%if (epsprint); print -depsc2 Figure3_intK5_intZ5.eps; end;
%if (epsprint); print -depsc2 Figure3_alt.eps; end;