% simulating OU process
clear all;

N = 1000000; 
dT = 0.25;
Zmean = 0;
mu = 0.25; % 1 - mu is persistence
sigma = 0.007;
% muini = gds; 
Zsim1 = zeros(N,1); 
Zsim2 = zeros(N,1); 
rng(100);  
shock = randn(N,1); 
shock(1,1) = 0;
mmu = -1 + mu;

% Zsim1 is used in fokker_planck.m
% Zsim2 is used in simulate.m (for Den Haan Errors)

for time = 1:N-1
    if time == 1
        Zsim1(time+1) = mu * dT * Zmean + (1 - mu * dT) * Zmean + sigma * shock(time) * sqrt(dT);
        Zsim2(time+1) = (1 - mmu * dT)^(-1) * (Zmean + sigma * shock(time) * sqrt(dT));
    else
        Zsim1(time+1) = mu * dT * Zmean + (1 - mu * dT) * Zsim1(time) + sigma * shock(time) * sqrt(dT);
        Zsim2(time+1) = (1 - mmu * dT)^(-1) * (Zsim2(time) + sigma * shock(time) * sqrt(dT));
    end
end

sigma/sqrt(1-(1-mu)^2) % theoretical moment?
std(Zsim1)
std(Zsim2)

figure;
plot(Zsim1);
hold on;
plot(Zsim2);