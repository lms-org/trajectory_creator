function [Jts] = otg_Jts(ps, T)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_Jts returns value of the functional Jt = int_0^T  (s'''(t))^2 dt for s
%   INPUT:
%       ps = coefficients of the polynomial (otg_ps) s.t.
%           s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
%       T = end time
%
%   OUTPUT:
%       Jts  = int_0^T  (s'''(t))^2 dt
%
% See also 

%precomputed formula
Jts = (192*T^3*ps(1)^2 + 144*T^2*ps(2)*ps(1) + 36*T*ps(2)^2);
end

