function [kappa] = kappaApprox(rbar,D,Regularize,kappaType,newtonApprox)
switch kappaType
    case 1
        kappa = max((rbar*D - rbar^3)/(1-rbar^2),Regularize);
        if newtonApprox
            kappa = kappaNewtonApprox(kappa, rbar,D, 2);
        end
    case 2
        kappa = max((D - 1) / 2*(1 - rbar), Regularize);
end
end

