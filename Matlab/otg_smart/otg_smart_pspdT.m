function [flag1, flag2, flag3, flagAll, ps, pd, T, TOL] = otg_smart_pspdT(absTOL, maxIter, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%This is the otg_smart subpackage which relies on the assumption that there
%is only one other vehicle on the road. This improves the optimation
%drastically in speed and accuracy. If this fails the otg_dumb subpackage
%should be used
%
%--------------------------------------------------------------------------
%
%otg_smart_opt_pspdT gives ps, pd, T for the optimal trajectory
%T
%   INPUT:
%
%       absTOL = (Tmax-Tmin) < absTOL --> return (the best T is guranteed
%           to lie within a ball with radius absTOL and middle Topt
%       maxIter = max. Number of Iterations (50 should be suffiecient)
%
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
%       dataVeh = (3x1) data of other vehivles on the road
%           dataVeh = [s0; v; I] s.t.
%               s0 = initial distance between obstacle car i and own car
%               v  = velocity of obstacle car i anlong the road
%               I  = -1/1: -1: on right lane, +1: on the left lane: this
%               must be the same lane as the start of the vehicle
%       safetyS = min. safety distance in s direction
%       safetyD = min. safety distance in D direction with sign:
%           if I == -1: d(t) must be bigger  than safetyD
%           if I ==  1: d(t) must be smaller than safetyD
%
%       kappa = curvature of the circle which describes the centerline
%       kappaMax = max. curvature of the road in xy coordinate system
%       aOrthMax = max. acceleration orthogonal to the trajectory
%
%   OUTPUT:
%       ps = coefficients of the polynomial in std. matlab notation i.e.
%           s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
%       pd = coefficients of the polynomial in std. matlab notation i.e.
%           d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
%           + pd(6)
%       T = end time
%       TOL = limit on the error in T
%
%       flag1 = flag for the subroutine otg_smart_objFun
%                1: all good
%               -1: the minimum velocity is negative!
%               -2: the safety point was not unique
%               -3: not considered case (by the programmer)
%       flag2 = flag for the subroutine otg_smart_opt_step
%                1: all good in this routine
%                0: no drivable trajectory
%               -1: driveability condition isn't as expected
%               -2: Ts not monotonically increasing
%       flag3 = flag for the otg_smart_optT subroutine
%                1: all good in this routine
%                0: no drivable trajectory
%               -1: number of iterations not sufficient
%                   -1.1: number of iterations not sufficient + sol. unique
%                   -1.2: number of iterations not sufficient + not unique
%               -2: no drivable solution was found with the given number of iterations
%               -10: something went wrong in one of the subroutines
%       flagAll = all together
%               1: all is good;
%               -1: something went somwhere wrong
%
%
% See also

%% init
ps = zeros(5,1);
pd = zeros(6,1);

flagAll = -1;

Tstart = dataVeh(1) / (abs((S(1)+S(3))/2 - 0.5*dataVeh(2)) + 0.001);

if (Tstart <= 1) || (Tstart >= 20)
    Tstart = 5; %just some hard-coded plausibility check
end

%% get the best T
[flag1, flag2, flag3, T, TOL] = otg_smart_optT(Tstart,absTOL, maxIter, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);

if (flag1 < 0) || (flag2 <= 0) || (flag3 <= 0)
    return; %if no return: all good
end

% get the coefficients
ps = otg_ps(S,T);
pd = otg_pd(D,T);

% see if not colliding
%too lazy, should be ok
[C, notD, coll, ~] = otg_smart_objFun(T, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);

if (notD(1) == 0) && (coll(1) == 0)
    flagAll = 1;
end







end