function [Jtd] = otg_Jtd(pd, T)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_Jtd returns value of the functional Jt = int_0^T  (d'''(t))^2 dt for d
%   INPUT:
%       pd = coefficients of the polynomial d (otg_pd) s.t.
%           d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
%           + pd(6)
%       T = end time
%
%   OUTPUT:
%       Jtd = int_0^T  (d'''(t))^2 dt
%
% See also 

%precomputed formula
Jtd = 36*T*pd(3)^2 + T^3*(192*pd(2)^2 + 240*pd(3)*pd(1)) + 720*T^5*pd(1)^2 + 144*T^2*pd(3)*pd(2) + 720*T^4*pd(2)*pd(1);

end

