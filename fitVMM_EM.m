function obj =  fitVMM_EM(X,Kmax,varargin)
[start,tol,maxIter,Regularize,debg,reps, kappaType] = process_options(varargin,'start',[],'tol',1e-12,'maxIter',300,'Regularize',1e-3,'debg',0,'reps',1, 'kappaType',1);
[N,D]= size(X);
%% Initialization using Spherical KMeans
% Inititalize mean, kappa and component mixing weight
if isempty(start)
    [Mu,Kappa,W] = InitParameters(X,Kmax,'kappaMin',Regularize,'debg',0,'reps', reps,'kappaType',kappaType);
else
    W = start.W;
    Mu = start.MU;
    Kappa = start.Kappa;
    Kmax = length(W);
end
% EM
iteration = 2;
loglike = -inf(maxIter,1);
converged = 0;
besselFunction = @logbesseliApprox1; %if change if numerical problems occur.

 while (~converged)
    % E-step
    [R,loglike(iteration)] = Expectation(X,W,Mu,Kappa, besselFunction);
    % M-step
    [W,Mu,Kappa] = Maximization(X,Kappa,R,Regularize,kappaType);
    % convergence check
    deltlike = loglike(iteration) - loglike(iteration-1);
    deltlike = abs(100*(deltlike/loglike(iteration-1)));
    if(deltlike < tol || iteration > maxIter)
        converged = 1;
        loglike(iteration+1:end) = [];
    end
    prt(debg,1,sprintf('########### EM Iteration: %d, LogLikelihood=%8.8f, Delta=',iteration,loglike(iteration)),deltlike);
    iteration = iteration + 1;
end
% Store results in object
obj.Iters = iteration-1;
obj.DistName = 'Mixture of vMF distributions';
obj.NDimensions = D;
obj.NSamples = N;
obj.NComponents = Kmax;
obj.PComponents = W;
obj.mu = Mu;
obj.Kappa = Kappa;
obj.E = R;
[~, idx] = max(R,[],2);
obj.Class =  idx;
obj.logL =  loglike(end);
end
