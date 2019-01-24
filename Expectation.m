function [R,llh] =  Expectation(X,W,Mu,Kappa)
% E-Step
[N, D]= size(X);
logNormalize  = log(W) + (D/2-1)*log(Kappa) - (D/2)*log(2*pi) - logbesseli(D/2-1,Kappa);
R  = X * (Mu'.*(ones(D,1)*Kappa));
R  = bsxfun(@plus,R,logNormalize);
T = logsumexp(R,2);
llh = sum(T)/N; % loglikelihood
R = exp(bsxfun(@minus,R,T));
end