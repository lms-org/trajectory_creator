function [flag1, flag2, Ts_new, Cs_new, notD_new, coll_new] = otg_smart_opt_step(Ts, Cs, notD, coll, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax)
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
%otg_smart_opt_step performs one step of the iterative search for the best
%T
%   INPUT:
%
%       Ts = [1x3] 3 time points of the current optim step
%       Cs = [1x3] 3 values of the cost function belonging to the Ts
%       notD = [1x3] 3 values of drivability cond. belonging to the Ts
%           = 0 <==> trajectory i is drivable
%           = 1 <==> trajectory i not drivable because of curvature
%           = 2 <==> trajectory i not drivable because of orth. acceleration
%           = 3 <==> trajectory i not drivable because of both above
%       coll = [1x3] 3 values of collision cond. belonging to the Ts
%           = 1 <==> trajectory i is colliding
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
%       Ts_new = [1x3] 3 time points of the next optim step
%       Cs_new = [1x3] 3 values of the cost function belonging to the Ts_new
%       notD_new = [1x3] 3 values of drivability cond. belonging to the
%               Ts_new
%           = 0 <==> trajectory i is drivable
%           = 1 <==> trajectory i not drivable because of curvature
%           = 2 <==> trajectory i not drivable because of orth. acceleration
%           = 3 <==> trajectory i not drivable because of both above
%       coll_new = [1x3] 3 values of collision cond. belonging to the
%               Ts_new
%           = 1 <==> trajectory i is colliding
%       flag1 = flag for the subroutine otg_smart_objFun 
%                1: all good
%               -1: the minimum velocity is negative!
%               -2: the safety point was not unique
%               -3: not considered case (by the programmer)
%       flag2 = flag for this routine
%                1: all good in this routine
%                0: no drivable trajectory
%               -1: driveability condition isn't as expected
%               -2: Ts not monotonically increasing 
%
% See also

%% init

flag2 = 1;
Ts_new = zeros(1,3);
Cs_new = zeros(1,3);
notD_new = zeros(1,3);
coll_new = zeros(1,3);
flag1 = 0;

%% check for order
if ~((Ts(1) < Ts(2))&&(Ts(2) < Ts(3)))
    flag2 = -2;
    return;
end

%% check if there is an error with the drivability
if (notD(3) > 0)
    %all points should be non drivable
    if (notD(1) == 0) || (notD(2) == 0)
        %something went wrong about the drivability
        flag2 = -1;
        return;
    end
end

%% check if there is an error with the collisions
if (coll(1) > 0)
    %all points should be colliding
    if (coll(2) == 0) || (coll(3) == 0)
        %something went wrong about the collsion
        flag2 = -3;
        return;
    end
end

%% check if no valid trajectory exists
if sum((notD > 0) & (coll > 0)) > 0
    %there is at least one traj wich is both not drivable and not colliding
    %--> no valid trajectory exists
    flag2 = 0;
    return
end


%% all non drivable
if (notD(1) > 0) && (notD(2) > 0) && (notD(3) > 0)
    Ts_new(1) = Ts(2);
    Ts_new(2) = Ts(3);
    Ts_new(3) = Ts(3) + (Ts(3)-Ts(1));
    
    [C3, notD3, coll3, flag1] = otg_smart_objFun(Ts_new(3), S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);
    
    if flag1 < 0
        %something in the subroutine went wrong
        return;
    end
    
    Cs_new(1) = Cs(2);
    Cs_new(2) = Cs(3);
    Cs_new(3) = C3;
    
    notD_new(1) = notD(2);
    notD_new(2) = notD(3);
    notD_new(3) = notD3;
    
    coll_new(1) = coll(2);
    coll_new(2) = coll(3);
    coll_new(3) = coll3;
    
    return;
end

%% all collidiong
if (coll(1) > 0) && (coll(2) > 0) && (coll(3) > 0)
    Ts_new(1) = Ts(1) - (Ts(3) - Ts(1)); %fibonacci growth
    Ts_new(2) = Ts(1);
    Ts_new(3) = Ts(2);
    
    [C1, notD1, coll1, flag1] = otg_smart_objFun(Ts_new(1), S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);
    
    if flag1 < 0
        %something in the subroutine went wrong
        return;
    end
    
    Cs_new(1) = C1;
    Cs_new(2) = Cs(1);
    Cs_new(3) = Cs(2);
    
    notD_new(1) = notD1;
    notD_new(2) = notD(1);
    notD_new(3) = notD(2);
    
    coll_new(1) = coll1;
    coll_new(2) = coll(1);
    coll_new(3) = coll(2);
    
    return;
end

%% Test all the cases
%workaround for cpp not offering a value inf: if one trajectory is not
%drivable or colliding set its c value to double the max

maxC = max(Cs);

Cs((notD>0) | (coll>0)) = maxC;

if (Cs(1) < Cs(2)) && (Cs(2) <= Cs(3))
    
    Ts_new(1) = Ts(1) - (Ts(3) - Ts(1)); %fibonacci growth
    Ts_new(2) = Ts(1);
    Ts_new(3) = Ts(2);
    
    [C1, notD1, coll1, flag1] = otg_smart_objFun(Ts_new(1), S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);
    
    if flag1 < 0
        %something in the subroutine went wrong
        return;
    end
    
    Cs_new(1) = C1;
    Cs_new(2) = Cs(1);
    Cs_new(3) = Cs(2);
    
    notD_new(1) = notD1;
    notD_new(2) = notD(1);
    notD_new(3) = notD(2);
    
    coll_new(1) = coll1;
    coll_new(2) = coll(1);
    coll_new(3) = coll(2);
    
    return;
end

if (Cs(1) >= Cs(2)) && (Cs(2) > Cs(3))
    Ts_new(1) = Ts(2);
    Ts_new(2) = Ts(3);
    Ts_new(3) = Ts(3) + (Ts(3)-Ts(1));%fibonacci growth
    
    [C3, notD3, coll3, flag1] = otg_smart_objFun(Ts_new(3), S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);
    
    if flag1 < 0
        %something in the subroutine went wrong
        return;
    end
    
    Cs_new(1) = Cs(2);
    Cs_new(2) = Cs(3);
    Cs_new(3) = C3;
    
    notD_new(1) = notD(2);
    notD_new(2) = notD(3);
    notD_new(3) = notD3;
    
    coll_new(1) = coll(2);
    coll_new(2) = coll(3);
    coll_new(3) = coll3;
    
    return;
end

if (Cs(1) >= Cs(2)) && (Cs(2) <= Cs(3))
    
    Ta = (Ts(1)+Ts(2))/2;
    Tb = (Ts(2)+Ts(3))/2;
    
    [Ca, notDa, colla, flag1] = otg_smart_objFun(Ta, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);
    if flag1 < 0
        %something in the subroutine went wrong
        return;
    end
    
    [Cb, notDb, collb, flag1] = otg_smart_objFun(Tb, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);
    if flag1 < 0
        %something in the subroutine went wrong
        return;
    end
    
    if (notDa > 0) || (colla > 0)
        Ca = maxC;
    end
    
    if (notDb > 0) || (collb > 0)
        Cb = maxC;
    end
    
    if (Cb < Cs(2))
        
        %out of the five points the minumum must be somewhere between the
        %three outter right ones
        
        Ts_new(1) = Ts(2);
        Ts_new(2) = Tb;
        Ts_new(3) = Ts(3);
        
        Cs_new(1) = Cs(2);
        Cs_new(2) = Cb;
        Cs_new(3) = Cs(3);
        
        notD_new(1) = notD(2);
        notD_new(2) = notDb;
        notD_new(3) = notD(3);
        
        coll_new(1) = coll(2);
        coll_new(2) = collb;
        coll_new(3) = coll(3);
        
        return;
    end
    
    if (Ca < Cs(2))
        
        %out of the five points the minumum must be somewhere between the
        %three outter left ones
        
        Ts_new(1) = Ts(1);
        Ts_new(2) = Ta;
        Ts_new(3) = Ts(2);
        
        Cs_new(1) = Cs(1);
        Cs_new(2) = Ca;
        Cs_new(3) = Cs(2);
        
        notD_new(1) = notD(1);
        notD_new(2) = notDa;
        notD_new(3) = notD(2);
        
        coll_new(1) = coll(1);
        coll_new(2) = colla;
        coll_new(3) = coll(2);
        
        return;
    end
    
    if (Ca >= Cs(2)) && (Cb >= Cs(2))
        
        %out of the five points the minumum must be somewhere between the
        %three middle ones
        
        Ts_new(1) = Ta;
        Ts_new(2) = Ts(2);
        Ts_new(3) = Tb;
        
        Cs_new(1) = Ca;
        Cs_new(2) = Cs(2);
        Cs_new(3) = Cb;
        
        notD_new(1) = notDa;
        notD_new(2) = notD(2);
        notD_new(3) = notDb;
        
        coll_new(1) = colla;
        coll_new(2) = coll(2);
        coll_new(3) = collb;
        
        return;
    end
    
end