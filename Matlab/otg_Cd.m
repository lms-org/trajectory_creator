function [Cd] = otg_Cd(pd, T, kj, kT)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_Cd returns value of the cost function Cd = kj*Jtd + kT*T for d
%   INPUT:
%       pd = coefficients of the polynomial d (otg_pd) s.t.
%           d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
%           + pd(6)
%       T = end time
%       kj = weight for the jerk functional
%       kT = weight for the time
%
%   OUTPUT:
%       Cd = kj*Jtd + kT*T
%
% See also 

%precomputed formula
Cd = kj*otg_Jtd(pd, T) + kT*T;
end




