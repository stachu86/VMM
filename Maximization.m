function [W,Mu,Kappa] = Maximization(X,Kappa,R,Regularize)
[N, K] = size(R);
[~, D] = size(X);
W = sum(R,1)./N;
Mu = R'*X;
for k=1:K
    normMu   = sqrt(Mu(k,:)*Mu(k,:)');
    rbar  = normMu/W(k);
    Mu(k,:)  = Mu(k,:)/normMu;
    Kappa(k) = max((rbar*D - rbar^3)/(1-rbar^2),Regularize);
end
end
