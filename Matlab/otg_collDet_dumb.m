function [collision] = otg_collDet_dumb(ps, pd, T, dataVeh, safetyS, safetyD, dt)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_collDet_dumb sees if a collision occurs (collision = 1)
%   INPUT:
%       ps = coefficients of the polynomial (otg_ps) s.t.
%           s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
%       pd = coefficients of the polynomial d (otg_pd) s.t.
%           d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
%           + pd(6)
%       T = end time
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
%   OUTPUT:
%       collsion = 0/1: 0: no collsion, 1: collsion
%
% See also


N = size(dataVeh, 2); %number of obstacle vehicles

tt = linspace(0, T, ceil(T/dt)); %samples
ss = polyval(ps, tt);
dd = polyval(pd, tt);

collision = 0;

for k = 1:length(tt)
    for i = 1:N %go over all vehicles
        
        Ii = dataVeh(3,i);%just for readability
        
        if (Ii == -1) && (dd(k) < safetyD)
            abstand = dataVeh(1,i) + dataVeh(2,i)*tt(k) - ss(k); %anfangsabstand + geschw. Auto * t - position unser Auto
            if abstand < safetyS
                collision = 1; %more than one collsion doesn't matter
                return;
            end
        end
        
        if (Ii ==  1) && (dd(k) > safetyD)
            abstand = dataVeh(1,i) + dataVeh(2,i)*tt(k) - ss(k); %anfangsabstand + geschw. Auto * t - position unser Auto
            if abstand < safetyS
                collision = 1; %more than one collsion doesn't matter
                return;
            end
        end
        
        
    end
end




