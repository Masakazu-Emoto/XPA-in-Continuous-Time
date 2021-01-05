clear all;

load ./XPA/Zsim_XPA.mat Zsim XPA_K Kds;
Zsim_XPA = Zsim;
Ksim_XPA = XPA_K;
Kss_XPA = Kds;

%load ./REITER/Zsim_REITER_rrhoTFP0.25.mat vAggregateTFP Kpath KSS;
load ./REITER/Zsim_REITER_rrhoTFP0.75.mat vAggregateTFP Kpath KSS;
Zsim_REITER = vAggregateTFP;
Ksim_REITER = Kpath;
Kss_REITER = KSS;

disp([Kss_XPA Kss_REITER]);

figure;
plot(Zsim_XPA);
hold on;
plot(Zsim_REITER);

figure;
plot(Ksim_XPA);
hold on;
plot(Ksim_REITER);