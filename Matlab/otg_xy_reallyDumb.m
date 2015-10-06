function [flag1, flag2, x, y] = otg_xy_reallyDumb(S, D, kj, kT, ks, kd, dT, Tmin, Tmax, dataVeh, safetyS, safetyD, dt, ds, m, kappa, y0, phi)
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
%               Ii  = -1/1: -1: on left lane, +1: on the right lane
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
%       y0 = distance from the center line along the y axis of the cosy of
%           the car
%       phi = angle between the x-axis of the car and the tangent on the
%           center line at the section between center line and y axis
%
%   OUTPUT:
%       flag1 = 1/0/-1
%           1:  everything ok
%           0:  solution was not unique. solution with minimal T is choosen
%           -1: no solution: the vectors are full of zeros
%       flag2 = 1/0/-1/-2
%           1:  everything ok
%           0.5: m too big, smaller one set
%           0:  solution was not unique. solution with minimal T is choosen
%           -1: the roots were imaginary
%           -2: no roots
%           -3: ds was too big
%       x = x in the car coordinate sytsem
%       y = y in the car coordinate sytsem
%
% See also

%% init
x = zeros(1,m);
y = zeros(1,m);

flag2 = 1;


%% coefficients
[flag1, ps, pd, T] = otg_pspdT_reallyDumb(S, D, kj, kT, ks, kd, dT, Tmin, Tmax, dataVeh, safetyS, safetyD, dt);

if flag1 == -1
    return;
end

if ((m-1)*ds > polyval(ps, T))
    flag2 = 0.5;
    m = floor(polyval(ps, T)/ds)-1;
    if m < 1
        flag2 = -3;
        return;
    end
end


%% calculate the times where s(t) = i*ds
tt = zeros(1,m);

for k = 2:m
    ps_new = ps - [0 0 0 0 (k-1)*ds];
    r = roots(ps_new);
    if length(r) == 0
        flag2 = -2;
        return;
    elseif length(r) > 1
        flag2 = 0;
        r = r(1);
    end
    
    if ~isreal(r)
        flag2 = -1;
        return;
    end
    tt(k) = r;
    
end

%% calc s and d at the times
ss = polyval(ps, tt);
dd = polyval(pd, tt);


XY = 1/kappa*[sin(kappa*ss);1-cos(kappa*ss)] + [-sin(kappa*ss);cos(kappa*ss)].*repmat(dd,2,1);

xy_car = zeros(2,m);

R = [cos(phi), -sin(phi); sin(phi) cos(phi)];


for k = 1:m
   xy_car(:,k) =  R*XY(:,k) + [-sin(phi);cos(phi)]*(cos(phi)*y0); 
end

x = xy_car(1,:);
y = xy_car(2,:);


end