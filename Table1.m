%Table 1
% set loadtemp = 0; (for baseline parameters) and savetime = 1; in main_XPA_v1.m, main_KS_v1.m,
% and main_REITER_v1.m before running this code
clear all;

cd ./REITER;

main_REITER_v1
load etime.mat etime1;
etimeall = zeros(10,1);
etimeall(1) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(2) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(3) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(4) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(5) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(6) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(7) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(8) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(9) = etime1;
save etimeall.mat etimeall;
main_REITER_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(10) = etime1;
save etimeall.mat etimeall;

cd ../

mean(etimeall)
std(etimeall)

cd ./XPA;

main_XPA_v1
load etime.mat etime1;
etimeall = zeros(10,1);
etimeall(1) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(2) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(3) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(4) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(5) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(6) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(7) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(8) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(9) = etime1;
save etimeall.mat etimeall;
main_XPA_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(10) = etime1;
save etimeall.mat etimeall;

cd ../

mean(etimeall)
std(etimeall)

cd ./KS;

main_KS_v1
load etime.mat etime1;
etimeall = zeros(10,1);
etimeall(1) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(2) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(3) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(4) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(5) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(6) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(7) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(8) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(9) = etime1;
save etimeall.mat etimeall;
main_KS_v1
load etime.mat etime1;
load etimeall.mat etimeall;
etimeall(10) = etime1;
save etimeall.mat etimeall;

cd ../

mean(etimeall)
std(etimeall)