function [kappa] = kappaNewtonApprox(kappa0, rbar,d, iter)
% kappaNewtonApprox this function performs few iterations of Newton's method
%   to make estimation of Kappa more accurate.
%
%   Based on "3.3 Truncated Newton Approximation" from
%   Sra, Suvrit. "A short note on parameter approximation for von
%   Mises-Fisher distributions: and a fast implementation of I s (x)." 
%   Computational Statistics 27.1 (2012): 177-190.
%
%   author: skacprza@agh.edu.pl

kappa = kappa0;
for i =1:iter
    A = besseli(d/2,kappa)./besseli(d/2-1,kappa);
    kappa = kappa - (A - rbar) ./ (1 - A.^2 - ((d-1)./kappa).*A);
end
if any(isnan(kappa)) || any(isinf(kappa)) 
    error('Numerical problem with Kappa!');
end
end

