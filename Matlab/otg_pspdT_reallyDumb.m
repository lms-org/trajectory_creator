function [flag, ps, pd, T] = otg_pspdT_reallyDumb(S, D, kj, kT, ks, kd, dT, Tmin, Tmax, dataVeh, safetyS, safetyD, dt)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_pspdT_reallyDumb returns the coefficients of the polynomial with the
%best cost function and no collision
%   INPUT:
%       
%       S = [v0, a0, v1] s.t.
%           s(0) = 0, s'(0) = v0, s''(0) = a0, s'(T) = v1
%       D = [d0, d0d, d0dd, d1] s.t.
%           d(0) = d0, d'(0) = d0d, d''(0) = d0dd, d'(T) = d1
%
%       kj = weight for the jerk functional
%       kT = weight for the time
%       ks = weight of the longitudinal weight function
%       kd = weight of the lateral weight function
%
%       dT = time interval between two possible end Times T
%       Tmin = minimal T
%       Tmax = maximal T
%
%       dataVeh = (3xN) data of other vehivles on the road
%           dataVeh(:,i) = [s0i; vi; Ii] s.t.
%               s0i = initial distance between obstacle car i and own car
%               vi  = velocity of obstacle car i anlong the road
%               Ii  = -1/1: -1: on right lane, +1: on the left lane
%       safetyS = min. safety distance in s direction
%       safetyD = min. safety distance in D direction with sign:
%           if Ii == -1: d(t) must be bigger  than safetyD
%           if Ii ==  1: d(t) must be smaller than safetyD
%       dt = sampling interval for the collision detection
%
%   OUTPUT:
%       flag = 1/0/-1
%           1:  everything ok
%           0:  solution was not unique. solution with minimal T is choosen
%           -1: no solution: the vectors are full of zeros
%       ps = coefficients of the polynomial in std. matlab notation i.e.
%           s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
%       pd = coefficients of the polynomial in std. matlab notation i.e.
%           d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
%           + pd(6)
%       T = end time
%
% See also

%% init for case of return
ps = zeros(1, 5);
pd = zeros(1, 6);
T = 0;
flag = 1;

%% init
n = ceil(Tmax-Tmin)/dT; %number samples in the T space

TT = linspace(Tmin, Tmax, n); %generate all Ts

PS = zeros(5, n); %all the coefficients pd

PD = zeros(6, n); %all the coefficients pd

C = zeros(1, n); %the cost function values
collisions = zeros(1, n); %collisions (0: no collision, 1: collision)

for k = 1:n
    PS(:, k) = otg_ps(S, TT(k)); %get coefficients
    PD(:,k)  = otg_pd(D, TT(k)); %get coefficients
    
    C(k) = otg_Ctot(PS(:, k), PD(:, k), TT(k), kj, kT, ks, kd); %compute total cost function
    
    collisions(k) = otg_collDet_dumb(PS(:, k), PD(:, k), TT(k),dataVeh, safetyS, safetyD, dt); %check for collision
end

ind_noCol = collisions == 0; %get all with no collisions

if sum(ind_noCol) == 0
    %no trajectory with no collision
    flag = -1;
    return;
end

ind_sol = find((C == min(C(ind_noCol))) & ind_noCol); %the second one could be useless but there is the slight chance that a colliding path has the same C value as the noncolliding minimal one

%% give solution
ps = PS(:,ind_sol(1));
pd = PD(:,ind_sol(1));
T = TT(ind_sol(1));

end