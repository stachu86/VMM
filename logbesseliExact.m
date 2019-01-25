function [logb] = logbesseliExact(nu,x)
exact = log(besseli(nu,x));
if any(isinf(exact)) || any(isnan(exact))
   error('Bessel function calculation problem!');  
end
logb = exact;