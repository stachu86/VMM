function obj =  fitVMM_CEM(X,Kmax,varargin)
%% Component Wise EM algorithm for Mixture of von Mises-Fisher distributions
%  "Unsupervised Learning of Mixture of von Mises-Fisher distributions
%  using Minimum Message Length"
%
%   obj = fitVMM_CEM(X,Kmax)
%   obj = fitVMM_CEM(X,Kmax,varargin)
%   optional inputs
%                 'Kmin' : the mininum number of compoent, default [1]
%                 'start' :  custom initialization start.W for mixing
%                 weights , start.MU for mean directions and start.Kappa for concentration parameters
%                 'Regularize' : the minimum Kappa allowed,
%                 'debg':  debug verbosity level , default [0] no info is printed
%                 'tol' :  covergence threshold [loglikelihood change between two consecutive iterations].
%
%   see also fitVMM_EM

%%
%#   $Author: Israel D. Gebru $    $Date: 2016/04/28 $    $Revision: 1.0 $
%#   $Modified: skacprza@agh.edu.pl $Date: 2019/02/06 $
%#   Copyright:
[Kmin,start,tol,Regularize,debg,reps,kappaType] = process_options(varargin,'Kmin',1,'start',[],'tol',1e-12,'Regularize',1e-3,'debg',0,'reps',10,'kappaType', 1);
[N,D]= size(X);

besselFunction = @logbesseliApprox1; %if change if numerical problems occur.

%% Initialization using Spherical KMeans
% Inititalize mean, kappa and component mixing weight
if isempty(start)
    [Mu,Kappa,W] = InitParameters(X,Kmax,'kappaMin',Regularize,'debg',debg, 'reps', reps, 'kappaType', kappaType);
else
    W = start.W;
    Mu = start.MU;
    Kappa = start.Kappa;
    Kmax = length(W);
end
%% MML
npars = D*(D + 1)/2; % number of free parameters to estimate
nparsover2 = npars / 2;
K = Kmax;
% Using the initial parameters, compute the log-likelihood
niter = 1;
loglike(niter) =  -inf; % stores the log-likelihood
dlength = -loglike(niter) + (nparsover2*sum(log(W))) + (nparsover2 + 0.5)*K*log(N); % description length
dl(niter) = dlength;  % stores the description length
K_history(niter) = K;  %  store the number of components
% the transitions vectors will store the iteration number at which components are killed.
% transitions1 stores the iterations at which components are  killed by the M-step,
% while transitions2 stores the iterations at which we force components to zero.
transitions1 = [];
transitions2 = [];
% minimum description length seen so far, and corresponding % parameter estimates
mindl = dl(niter);
bestW = W;
bestMu = Mu;
bestKappa = Kappa;
bestK = K;
k_cont = 1;    % auxiliary variable for the outer loop
%%
while(k_cont)  % the outer loop will take us down from Kmax to Kmin components
    cont=1;        % auxiliary variable of the inner loop
    while(cont)    % this inner loop is the component-wise EM algorithm with the
        prt(debg,3,sprintf('k = %2d,  minestpp = %0.5g @iter=', K, min(W)),niter);
        % we begin at component 1
        comp = 1;
        while comp <= K
            [R,~] = Expectation(X,W,Mu,Kappa,besselFunction);
            % now we perform the standard M-step for Mean and Kappa
            [W,Mu,Kappa] = Maximization(X,Kappa,R,Regularize,kappaType);
            % this is the special part of the M step that is able to kill components
            W(comp) = max(sum(R(:,comp))-nparsover2,0)/N;
            W = W/sum(W);
            % this is an auxiliary variable that will be used to signal the killing of the current component being updated
            killed = 0;
            % do some book-keeping if the current component was killed
            if W(comp)==0
                prt(debg,2,'component killed..',comp);
                killed = 1;
                % register that at the current iteration a component was killed
                transitions1 = [transitions1,niter];                
                Kappa(comp) = [];
                Mu(comp,:) = [];
                W(comp) = [];                
                % since we killed a component, k must decrease
                K = K-1;
            end % end of W(comp)==0
            % if the component was not killed
            if killed==0
                comp = comp + 1;
            end
            % if killed==1, it means the in the position "comp", we now have a component that was not yet visited in this sweep,
            % and so all we have to do is go back to the M setp without
            % increasing "comp". This is  doing start EM steps
        end % this is the end of the innermost "while comp <= k" loop which cycles through the components
        % increment the iterations counter
        niter = niter + 1;
        [~,llh] = Expectation(X,W,Mu,Kappa,besselFunction); %perform E-step
        loglike(niter) = llh;
        % compute and store the description length and the current number of components
        dlength = -loglike(niter) + (nparsover2*sum(log(W))) +  (nparsover2 + 0.5)*K*log(N);
        dl(niter) = dlength;
        K_history(niter) = K;
        % compute the change in loglikelihood to check if we should stop
        deltlike = loglike(niter) - loglike(niter-1);
        prt(debg,1,sprintf('########### iter: %d, deltaloglike = %0.12f%% ,K=%d, BestK=',niter,abs(100*(deltlike/loglike(niter-1))),K),bestK);
        if abs(100*(deltlike/loglike(niter-1))) < tol
            % if the relative change in loglikelihood is below the tolerance threshold, we stop
            cont=0;
        end
    end % this end is of the inner loop: "while(cont)"
    % now check if the latest description length is the best if it is, we store its value and the corresponding estimates
    if dl(niter) < mindl
        bestW = W;
        bestMu = Mu;
        bestKappa = Kappa;
        bestK = K;
        mindl = dl(niter);
    end
    
    % at this point, we may try smaller mixtures by killing the component with the smallest mixing probability
    % and then restarting CEM2 as long as K is not yet at Kmin
    if K>Kmin
        [~, indminw] = min(W);
        Kappa(indminw) = [];
        Mu(indminw,:) = [];
        W(indminw) = [];
        K = K-1;
        % we renormalize the mixing probabilities after killing the component
        W = W./sum(W);
        % and register the fact that we have forced one component to zero
        transitions2 = [transitions2, niter];
        % increment the iterations counter
        niter = niter + 1;
        % compute the loglikelihhod function and the description length
        [~,llh] = Expectation(X,W,Mu,Kappa,besselFunction); %perform E-step
        loglike(niter) = llh;
        dlength = -loglike(niter) + (nparsover2*sum(log(W))) +  (nparsover2 + 0.5)*K*log(N);
        dl(niter) = dlength;
        K_history(niter) = K;
    else
        %if k is not larger than kmin, we must stop
        k_cont = 0;
    end
end % this is the end of the outer loop "while(k_cont)"
%%
% to merge exactly similar component. This happen when we have points with strong concentrated in one area
% the algorithm will fail to annihilate one of them, thus we do it here
mu = bestMu;
D = inf(bestK,bestK);
D_zero = zeros(bestK,bestK);
for i=1:bestK
    for j=i+1:bestK
        D(i,j) = roundn(pdist2(mu(i,:),mu(j,:)),-3);
    end
end
[r,~,~] = find(D==D_zero);
%remove one of them
if ~isempty(r)
    prt(debg,2,'Identical comp, removing one of them!',1);
    rr = unique(r);
    bestK = bestK-length(rr);
    bestW(rr)=[];
    bestMu(rr,:)=[];
    bestKappa(rr)=[];
end
% need to recompute
[R,llh] = Expectation(X,bestW,bestMu,bestKappa,besselFunction);
mindl = -llh + (nparsover2*sum(log(bestW))) +  (nparsover2 + 0.5)*bestK*log(N);
prt(debg,1,sprintf('########### Converged!!!,MML=%6.6f, BestK=',mindl),bestK);
% Store results in object
obj.Iters = niter;
obj.DistName = 'Mixture of vMF distribution';
obj.NDimensions = D;
obj.NSamples = N;
obj.NComponents = bestK;
obj.PComponents = bestW;
obj.mu = bestMu;
obj.Kappa = bestKappa;
obj.E = R;
[~, idx] = max(R,[],2);
obj.Class =  idx;
obj.logL =  loglike(niter);
obj.dl = dl;
obj.mindl = mindl;
end
