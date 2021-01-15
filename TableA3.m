% Table A4
clear all;
format short;

kub = 0.2;
klb = 0.2;
intK = 3;
intZ = 3;
vsig = [0.0001 0.001 0.007 0.01 0.03 0.05]';

DH_KS = zeros(size(vsig,1),3);
DH_KS_FVHN = zeros(size(vsig,1),3);  % UpwindKZ=0, KFEnoKZ=0
DH_KS_FVHN1 = zeros(size(vsig,1),3); % UpwindKZ=0, KFEnoKZ=1
DH_KS_FVHN2 = zeros(size(vsig,1),3); % UpwindKZ=1, KFEnoKZ=0
for i = 1:size(vsig,1)
    
    sigma = vsig(i);
    eval(sprintf('load ./KS/CT_KS_sigma%1.4f_kub%1.2f_klb%1.2f_intK%d_intZ%d.mat',sigma,kub,klb,intK,intZ));
    DH_KS(i,1) = DH_Max;
    DH_KS(i,2) = DH_Mean;
    DH_KS(i,3) = Y_R2;
    eval(sprintf('load ./KS/CT_KS_sigma%1.4f_kub%1.2f_klb%1.2f_intK%d_intZ%d_FVHN.mat',sigma,kub,klb,intK,intZ));
    DH_KS_FVHN(i,1) = DH_Max;
    DH_KS_FVHN(i,2) = DH_Mean;
    DH_KS_FVHN(i,3) = Y_R2;
    eval(sprintf('load ./KS/CT_KS_sigma%1.4f_kub%1.2f_klb%1.2f_intK%d_intZ%d_FVHN1.mat',sigma,kub,klb,intK,intZ));
    DH_KS_FVHN1(i,1) = DH_Max;
    DH_KS_FVHN1(i,2) = DH_Mean;
    DH_KS_FVHN1(i,3) = Y_R2;    
    eval(sprintf('load ./KS/CT_KS_sigma%1.4f_kub%1.2f_klb%1.2f_intK%d_intZ%d_FVHN2.mat',sigma,kub,klb,intK,intZ));
    DH_KS_FVHN2(i,1) = DH_Max;
    DH_KS_FVHN2(i,2) = DH_Mean;
    DH_KS_FVHN2(i,3) = Y_R2;

end

round(DH_KS',3)
round(DH_KS_FVHN',3)
round(DH_KS_FVHN1',3)
round(DH_KS_FVHN2',3)