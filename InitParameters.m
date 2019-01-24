function [mu,kappa,w] = InitParameters(X,K,varargin)
%% InitParameters
% Initialize Mean Vector and Concentration parameters of Mixture of von
% Mises-Fisher distributions. First the points in the N-by-D data matrix X are partitioned into K clusters
% using Spherical Kmeans.
% 
%   [mu,kappa,w] = InitParameters(X, K) returns the cluster centroid locations in the mu (K-by-D matrix center.
%   [mu,kappa,w] =  InitParameters(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the iterative algorithm used by KMEANS.  
%   Parameters are:
%
%   'kappaMin'  -  the minimum concentration parameter; Default is 1000.
%   'reps' - Number of times to repeat the clustering, each with a new set of initial centroids. Default is 10. 
%   'debg' = Verbosity level.  0 [default]

%#   $Author: Israel D. Gebru $    $Date: 2016/04/28 $    $Revision: 1.0 $
%#   Copyright:

%% Initialize parameters using Spherical Kmeans 
[kappaMin, debg,reps] = process_options(varargin,'kappaMin',100,'debg',0,'reps',10);
[n,d] = size(X);
%%
bestlabels = kmeans(X,K,'Distance','cosine','Replicate',reps);

% initializing kappa, mixing weight 
kappa = kappaMin*ones(1,K);
w = ones(1,K).*(1/n);
mu = zeros(K,size(X,2));
for k = 1:K
    idx = bestlabels == k;
    if sum(idx)>0
        w(k) = sum(idx)/n;   
        mu(k,:) = sum(X(idx,:));
        normMu   = norm(mu(k,:));
        rbar  = normMu/(n*w(k));
        mu(k,:)  = mu(k,:)./normMu;
        kappa(k) = max((rbar*d - rbar^3)/(1-rbar^2), kappaMin);
    else
        prt(debg,3,'Empty Cluster found!',[]);
    end
end
end