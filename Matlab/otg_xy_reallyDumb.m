function [flag, x, y, T] = otg_xy_reallyDumb(v1, d1, kj, kT, ks, kd, dT, Tmin, Tmax, dataVeh, safetyS, safetyD, dt, m, kappa, y0, phi, vx0, ax0, w)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_pspdT_reallyDumb returns the points specified by x, y on the
%trajectory in the car coordinate system
%   INPUT:
%
%       v1 = velocity in s direction at the end of the trajectory
%       d1 = distance from the center line at the end of the trajectory
%
%       kj = weight for the jerk functional
%       kT = weight for the time
%       ks = weight of the longitudinal weight function
%       kd = weight of the lateral weight function
%
%       dT = time interval between two possible end Times T
%       Tmin = minimal T > 0
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
%       ds = distance in s direction between two points
%       m = number of points wanted
%
%       kappa = curvature of the center line
%       phi = angle between the x-axis of the car and the tangent on the
%           center line at the section between center line and y axis
%
%   OUTPUT:
%       flag = 1/0/-1
%           1:  everything ok
%           0:  solution was not unique. solution with minimal T is choosen
%           -1: no solution: the vectors are full of zeros
%       x = x in the car coordinate sytsem
%       y = y in the car coordinate sytsem
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
[flag, ps, pd, T] = otg_pspdT_reallyDumb(S, D, kj, kT, ks, kd, dT, Tmin, Tmax, dataVeh, safetyS, safetyD, dt);

if flag == -1
    return;
end


%% calculate tt
tt = linspace(0,T,m);

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