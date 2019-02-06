function [logb] = logbesseli(nu,x)
% log of the Bessel function, extended for large nu and x
% approximation from Eqn 9.7.7 of Abramowitz and Stegun
% http://www.math.sfu.ca/~cbm/aands/page_378.htm
frac = x/nu;
square = 1 + frac.^2;
root = sqrt(square);
eta = root + log(frac) - log(1+root);
approx = - log(sqrt(2*pi*nu)) + nu*eta - 0.25*log(square);
logb = approx;

% alternative less accurate approximation from Eqn 9.6.7
% this gives better results on the Classic400 dataset!
% logb = nu*log(x/2) - gammaln(nu+1);





end