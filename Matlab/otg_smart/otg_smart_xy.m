function [flag1, flag2, flag3, flagAll, x, y, T, TOL] = otg_smart_xy(absTOL, maxIter, v1, d1, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappaMax, aOrthMax, m, kappa, y0, phi, vx0, ax0, w)
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
%otg_smart_opt_xy gives the coordinates x and y of the optimal trjectory in
%the car coordinate system
%   INPUT:
%
%       absTOL = (Tmax-Tmin) < absTOL --> return (the best T is guranteed
%           to lie within a ball with radius absTOL and middle Topt
%       maxIter = max. Number of Iterations (50 should be suffiecient)
%
%
%       v1 = velocity in s direction at the end of the trajectory
%       d1 = distance from the center line at the end of the trajectory
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
%       kappaMax = max. curvature of the road in xy coordinate system
%       aOrthMax = max. acceleration orthogonal to the trajectory
%
%       m = number of points wanted
%
%       kappa = curvature of the center line
%       phi = angle between the x-axis of the car and the tangent on the
%           center line at the section between center line and y axis
%
%   OUTPUT:
%       x = x in the car coordinate sytsem
%       y = y in the car coordinate sytsem
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
x = zeros(1,m);
y = zeros(1,m);

%% converte the state to S, D
S = [cos(phi)*vx0,-sin(phi)*w*vx0+cos(phi)*ax0, v1]; %really dirty

cosphi_d = -sin(phi)*w;
cosphi_dd = -w^2*cos(phi);
y0_d = -sin(phi)*vx0;
y0_dd = -(cos(phi)*w*vx0+sin(phi)*ax0);

D = [-cos(phi)*y0,-(cosphi_d*y0+cos(phi)*y0_d),-(cosphi_dd*y0+2*cosphi_d*y0_d+cos(phi)*y0_dd),d1];


%% coefficients
[flag1, flag2, flag3, flagAll, ps, pd, T, TOL] = otg_smart_pspdT(absTOL, maxIter, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax);

if flagAll == -1
    return;
end


%% calculate tt
tt = linspace(0,T(1),m);

%% calc s and d at the times
ss = polyval(ps, tt);
dd = polyval(pd, tt);

if kappa >0
    XY = 1/kappa*[sin(kappa*ss);1-cos(kappa*ss)] + [-sin(kappa*ss);cos(kappa*ss)].*repmat(dd,2,1);
else
    XY = [ss;zeros(1, m)] + repmat([0;1],1,m).*repmat(dd,2,1);
end

xy_car = zeros(2,m);

R = [cos(phi), -sin(phi); sin(phi) cos(phi)];


for k = 1:m
    xy_car(:,k) =  R*XY(:,k) + [-sin(phi);cos(phi)]*(-D(1));
end

x = xy_car(1,:);
y = xy_car(2,:);





end