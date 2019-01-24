function [W,Mu,Kappa] = Maximization(X,Kappa,R,Regularize)
[n, K] = size(R);
[~, D] = size(X);
W = sum(R,1)./n;
Mu = R'*X;
for k=1:K
    normMu   = norm(Mu(k,:));
    rbar  = normMu/(n*W(k));
    Mu(k,:)  = Mu(k,:)/normMu;
    Kappa(k) = max((rbar*D - rbar^3)/(1-rbar^2),Regularize);
end
end
