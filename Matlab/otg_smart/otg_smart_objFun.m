function [Ctot, notD, coll, flag] = otg_smart_objFun(T, S, D, kj, kT, ks, kd, dataVeh, safetyS, safetyD, kappa, kappaMax, aOrthMax)
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
%otg_smart_objFun returns for given T and parameters the value of the cost
%function C_tot and indicators if the trajectory corresponding to T is
%drivable
%   INPUT:
%
%       T = [1xn] for n different end times T(i)
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
%       C_tot = [1xn]: cost function values for the different T(i)
%       notD = [1xn] = 0/1/2/3 indicator if the trajectory is not drivable: notD(i)
%           = 0 <==> trajectory i is drivable
%           = 1 <==> trajectory i not drivable because of curvature
%           = 2 <==> trajectory i not drivable because of orth. acceleration
%           = 3 <==> trajectory i not drivable because of both above
%       coll = [1xn]: indicator if the trajectory is colliding: coll(i)
%           = 1 <==> trajectory i is colliding
%       flag
%           1: all good
%           -1: the minimum velocity is negative!
%           -2: the safety point was not unique
%           -3: not considered case (by the programmer)
%
% See also

%% init
n = length(T);
flag = 1;

Ctot = zeros(1,n);
notD = zeros(1,n);
coll = zeros(1,n);


%% Get the values of ps, pd for the Ts
PS = zeros(5,n);
PD = zeros(6,n);

for i = 1:n
    PS(:,i) = otg_ps(S,T(i));
    PD(:,i) = otg_pd(D,T(i));
end

%% Calculate the costs

for i = 1:n
    Ctot(i) = otg_Ctot(PS(:,i),PD(:,i),T(i), kj, kT, ks, kd);
end


%% See if drivable: We use a very rough approximation here!!

for i = 1:n
    %LENKEINSCHLAG
    d_d = polyder(PD(:,i)); %first derivative w.r.t. t of d(t)
    d_dd = polyder(d_d); %second derivative w.r.t. t of d(t)
    d_ddd = polyder(d_dd); %third derivative to find the max of the second derivative
    r = roots(d_ddd)';%get all roots of THIRD derivative
    
    ind = imag(r) == 0;%only real roots
    r = r(ind);
    
    ind = (r>0) & (r < T(i)); %only roots in interval
    r = r(ind);
    
    ts_d = [0 r T(i)]; %also max at the border
    d_dd_ts = polyval(d_dd, ts_d); %HERE ONE NEEDS THE SECOND DERIVATIVE
    
    d_dd_max = max(d_dd_ts); %get the maximal value of the polynomial on the interval
    
    
    s_d = polyder(PS(:,i)); %first derivative w.r.t. t of s(t) = v(t)
    s_dd = polyder(s_d); %second derivative to find the
    r = roots(s_dd)';%get all roots of first derivative
    
    ind = imag(r) == 0;%only real roots
    r = r(ind);
    
    ind = (r>0) & (r < T(i)); %only roots in interval
    r = r(ind);
    
    ts_s = [0 r T(i)]; %also min at the border
    s_d_ts = polyval(s_d, ts_s);
    
    s_d_min = min(s_d_ts); %get the min value of the polynomial on the interval
    s_d_max = max(s_d_ts); %get the maximal value of the velocity in s direction
    
    if (s_d_min) <= 0
        %something went terrible wrong. we are driving backwards
        flag = -1;
        return;
    end
    
    %now the bound for the curvature can be calculated
    kappa_xy_max = abs(kappa) + d_dd_max/s_d_min;
    
    if kappa_xy_max > kappaMax
        %the curvature could be (according to our approximation) to big
        notD(i) = 1;
    end
    
    
    %MAXIMALE ORTH: BESCHLEUNIGUNG
    r = roots(d_dd)';%get all roots of SECOND derivative
    
    ind = imag(r) == 0;%only real roots
    r = r(ind);
    
    ind = (r>0) & (r < T(i)); %only roots in interval
    r = r(ind);
    
    ts_d_2 = [0 r T(i)]; %also max at the border
    d_d_ts = polyval(d_d, ts_d_2); %HERE ONE NEEDS THE SECOND DERIVATIVE
    
    d_d_max = max(d_d_ts); %get the maximal value of the polynomial on the interval
    
    vx_max = sqrt(s_d_max^2+d_d_max^2);
    
    if vx_max^2*kappa_xy_max > aOrthMax
        %the orthogonal acceleration is too much for the grip
        notD(i) = notD(i)+2; %must be 2. 2+1 = 3 (all), 0 + 2= 2 (only acc.)
    end
end



%% check for collision
for i = 1:n
    %s(t) should be an increasing function. This is guranteed because we
    %calculate the min of the 1st derivative above and throw an error if
    %this is non positive
    
    abstandMinusSafety_poly = [0 0 0 dataVeh(2) dataVeh(1)-safetyS]' - PS(:,i);
    
    r = roots(abstandMinusSafety_poly);
    
    ind = imag(r) == 0;%only real roots
    r = r(ind);
    
    ind = (r>0) & (r < T(i)); %only roots in interval
    r = r(ind);
    
    if (isempty(r)) && (polyval(abstandMinusSafety_poly, T(i)/2) > 0)
        %all good no collsion possible
    elseif (isempty(r)) && (polyval(abstandMinusSafety_poly, T(i)/2) <= 0)
        %fix collsion
        coll(i) = 1;
    elseif (length(r) == 1) && (polyval(abstandMinusSafety_poly, r(1)/2) > 0)
        %one unique possible collision possible
        d_r = polyval(PD(:,i),r(1));
        
        I = dataVeh(3);
        
        if (I == -1) && (d_r < safetyD)
            coll(i) = 1;
        end
        
        if (I ==  1) && (d_r > safetyD)
            coll(i) = 1;
        end
    elseif (length(r) == 1) && (polyval(abstandMinusSafety_poly, r(1)/2) > 0)
        %there is a collision directly at the start
        coll(i) = -1;
    else
        flag = -3;
        return;
    end
end