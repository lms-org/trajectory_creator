%otg_smart test script
clear all;
close all;

%% test of the objective Function
dataVeh = [7; 0; -1];

D = [-0.2, 0, 0, 0.2];
S = [1, 0, 5];

kappa = 0.2;

T = linspace(1.5, 15, 500);

kj = 10;
kT = 1;
ks = 1;
kd = 1;

safetyS = 1.5;
safetyD = 0.1;

kappaMax = 1;

aOrthMax = 10;

[Ctot, notD, coll, flag] = otg_smart_objFun(T, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);

flag

subplot(3,1,1); plot(T,Ctot); subplot(3,1,2); stairs(T,notD); subplot(3,1,3); stairs(T,coll)
