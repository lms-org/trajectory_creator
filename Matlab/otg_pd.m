function [pd] = otg_pd(D, T)
%This is the OTG = Optimal Trajectory Generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_pd returns the coefficients of the unique 5th order Polynomial for d(t)
%   INPUT:
%       D = [d0, d0d, d0dd, d1] s.t.
%           d(0) = d0, d'(0) = d0d, d''(0) = d0dd, d'(T) = d1
%       T = end time
%
%   OUTPUT:
%       pd = coefficients of the polynomial in std. matlab notation i.e.
%           d(t) = pd(1)*t^5 + pd(2)*t^4 + pd(3)*t^3 + pd(4)*t^2 + pd(5)*t
%           + pd(6)
%
% See also 

%init
pd = zeros(1, 6);

%precomputed formulas
pd(1) = -(D(3)*T^2 + 6*D(2)*T + 12*D(1) - 12*D(4))/(2*T^5);
pd(2) = (3*D(3)*T^2 + 16*D(2)*T + 30*D(1) - 30*D(4))/(2*T^4);
pd(3) = -(3*D(3)*T^2 + 12*D(2)*T + 20*D(1) - 20*D(4))/(2*T^3);
pd(4) = D(3)/2;
pd(5) = D(2);
pd(6) = D(1);



end

