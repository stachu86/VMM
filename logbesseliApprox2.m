function [logb] = logbesseliApprox2(nu,x)
% alternative less accurate approximation from Eqn 9.6.7
% this gives better results on the Classic400 dataset!
logb = nu*log(x/2) - gammaln(nu+1);
% [bessel,flags] = besseli(nu,x);
% if any(flags ~= 0) || any(bessel == Inf)
%     besselproblem = [x, bessel, flags];
% end
% logb = bessel;
% nz = find(bessel > 0);
% z = find(bessel == 0);
% logb(nz) = log(bessel(nz));
% logb(z) = nu*log(x(z)/2) - gammaln(nu+1);