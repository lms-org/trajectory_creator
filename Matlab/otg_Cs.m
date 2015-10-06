function [Cs] = otg_Cs(ps, T, kj, kT)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_Cs returns value of the cost function Cs = kj*Jts + kT*T for s
%   INPUT:
%       ps = coefficients of the polynomial (otg_ps) s.t.
%           s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
%       T = end time
%       kj = weight for the jerk functional
%       kT = weight for the time
%
%   OUTPUT:
%       Cs = kj*Jts + kT*T
%
% See also 

%precomputed formula
Cs = kj*otg_Jts(ps, T) + kT*T;

end

