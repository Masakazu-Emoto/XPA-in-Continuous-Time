% Table A5
clear all;
format short;

mu = 0.75
kub = 0.2;
klb = 0.2;
intK = 3;
intZ = 3;
vsig = [0.0001 0.001 0.007 0.01 0.03 0.05]';

DH_KS = zeros(size(vsig,1),2);
DH_XPA = zeros(size(vsig,1),2);
DH_REITER = zeros(size(vsig,1),2);
for i = 1:size(vsig,1)
    
    sigma = vsig(i);
    eval(sprintf('load ./KS/CT_KS_mu%1.2f_sigma%1.4f_kub%1.2f_klb%1.2f_intK%d_intZ%d.mat',mu,sigma,kub,klb,intK,intZ));
    DH_KS(i,1) = DH_Max;
    DH_KS(i,2) = DH_Mean;
    eval(sprintf('load ./XPA/CT_XPA_mu%1.2f_sigma%1.4f_kub%1.2f_klb%1.2f_intK%d_intZ%d.mat',mu,sigma,kub,klb,intK,intZ));
    DH_XPA(i,1) = DH_Max;
    DH_XPA(i,2) = DH_Mean;
    eval(sprintf('load ./REITER/CT_REITER_mu%1.2f_sigma%1.4f.mat DH_Max DH_Mean',mu,sigma));
    DH_REITER(i,1) = DH_Max;
    DH_REITER(i,2) = DH_Mean;
    
end

% round(DH_XPA',3)
% round(DH_KS',3)
% round(DH_REITER',3)
[DH_XPA';
DH_KS';
DH_REITER']