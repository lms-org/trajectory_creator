function [flag1, flag2, flag3, Topt, TOL] = otg_smart_optT(Tstart,absTOL, maxIter, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax)
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
%otg_smart_optT performs the optimization in T
%T
%   INPUT:
%
%       Tstart = start value for T
%       absTOL = (Tmax-Tmin) < absTOL --> return (the best T is guranteed
%           to lie within a ball with radius absTOL and middle Topt
%       maxIter = max. Number of Iterations (50 should be suffiecient)
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
%       Topt = optimal time for the trajectory
%       TOL = limit on the error in T
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
%       flag3 = flag for this routine
%                1: all good in this routine
%                0: no drivable trajectory
%               -1: number of iterations not sufficient
%                   -1.1: number of iterations not sufficient + sol. unique
%                   -1.2: number of iterations not sufficient + not unique
%               -2: no drivable solution was found with the given number of iterations
%               -10: something went wrong in one of the subroutines
%
%
% See also

%% init
TOL = -100;
Topt = 0;
flag2 = -100;
flag3 = 1;

Ts = [0.5, 1, 1.5] * Tstart;
[Cs, notD, coll, flag1] = otg_smart_objFun(Ts, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);

for iter = 1:maxIter
    % call the subroutine for one step
    [flag1, flag2, Ts, Cs, notD, coll] = otg_smart_opt_step(Ts, Cs, notD, coll, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);
    
    if (flag1 < 0) || (flag2 <= 0)
        flag3 = -10;
        return;
    end
    
    if abs(Ts(3)-Ts(1)) < absTOL && (notD(2) == 0) && (coll(2) == 0)
        flag3 = 1;
        Topt = Ts(2);
        TOL = Ts(3) - Ts(1);
        return;
    end
    
    
    
end

%% number of iterations was not sufficient
indOk = (notD == 0) & (coll == 0);

if sum(indOk) > 0
    %at least one drivable path
    if sum(indOk) == 1
        %the drivable one is unique
        flag3 = -1.1;
        Topt = Ts(indOk);
        TOL = Ts(3) - Ts(1);
        return
    end
    
    % the drivable path is not unique, just choose the first one
    flag3 = -1.1;
    indOk = find(indOk);
    Topt = Ts(indOk(1));
    TOL = Ts(3) - Ts(1);
    return;
end

%no drivable solution was found with the given number of iterations
flag3 = -2;
return


end