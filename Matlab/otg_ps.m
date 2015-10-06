function [ps] = otg_ps(S, T)
%This is the OTG = Optimal trajectory generation package. This is based on
%the paper "Optimal Trajectory Generation for Dynamic Street Scenarios in a
%Frenet Frame"(2010) by Werling et. al.
%
%--------------------------------------------------------------------------
%
%otg_ps returns the coefficients of the unique 4th order Polynomial for s(t)
%   INPUT:
%       S = [v0, a0, v1] s.t.
%           s(0) = 0, s'(0) = v0, s''(0) = a0, s'(T) = v1
%       T = end time
%
%   OUTPUT:
%       ps = coefficients of the polynomial in std. matlab notation i.e.
%           s(t) = ps(1)*t^4 + ps(2)*t^3 + ps(3)*t^2 + ps(4)*t + ps(5)
%
% See also 

%init
ps = zeros(1, 5);


%formulas precomputed
ps(1) = 1.0./T.^3.*(S(1).*2.0-S(3).*2.0+T.*S(2)).*(1.0./4.0);
ps(2) = 1.0./T.^2.*(S(1).*3.0-S(3).*3.0+T.*S(2).*2.0).*(-1.0./3.0);
ps(3) = S(2).*(1.0./2.0);
ps(4) = S(1);
ps(5) = 0;


end

