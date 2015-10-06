function [Ctot] = otg_Ctot(ps, pd, T, kj, kT, ks, kd)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_Ctot returns value of the cost total function Ctot = ks*Cs + kd*Cd
%   INPUT:
%       ps = coefficients of the polynomial (otg_ps) s.t.
%           s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
%       pd = coefficients of the polynomial d (otg_pd) s.t.
%           d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
%           + pd(6)
%       T = end time
%       kj = weight for the jerk functional
%       kT = weight for the time
%       ks = weight of the longitudinal weight function
%       kd = weight of the lateral weight function
%
%   OUTPUT:
%       Ctot = ks*Cs + kd*Cd
%
% See also 

Ctot = ks*otg_Cs(ps, T, kj, kT) + kd*otg_Cd(pd, T, kj, kT);
end




