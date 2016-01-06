function [GMM,loglikelihood] = OnlineGMMEM(GMM0,dataSource,nIter,miniBatchSize,outputFileName,T0,alpha,firstBatchSize,removeFlatPatches)
% Learn a GMM model using an online version of the EM algorithm. 
% Parameters: 
%   GMM0: An initlization GMM or the number of components you want to learn
%   dataSource: the function you use to extract samples from the images
%   nIter: number of EM iterations
%   miniBatchSize: number of patches you use to estimate GMM after 10 iterations
%   outputFileName: name of output mat file containing GMM parameters
%   T0:
%   alpha: 
%   firstBatchSize: number of patches you use to estimate GMM for first 10 iterations
%   removeFlatPatches: true if you don't want to include flat patches

% sample call: NewGMM = OnlineGMMEM(100,@(N) removeDC(Rand3DPatchesFromImagesCell(N,PatchSize,Images)),10000,1500000,filename,2,0.9,1500000,true);

VERBOSE = true; 

%% =================== INITIALIZE PARAMETERS =================== %%

% if you pass in a struct with the GMM parameters, initilize it to this
if isstruct(GMM0)
    GMM = GMM0;
    K = GMM.nmodels;
    
% otherwise initilize GMM struct to have GMM0 components  
else
    GMM = [];
    K = GMM0;
    GMM.nmodels = K;
    GMM.mixweights = zeros(1,K);
    GMM.covs = zeros(size(dataSource(1),1),size(dataSource(1),1),K);
end

loglikelihood = zeros(1,nIter);

% assign parameter values that were not passed in 
if ~exist('T0','var')
    T0 = 500;
end
if ~exist('alpha','var')
    alpha = 0.6;
end
if ~exist('FirstBatchSize','var')
    firstBatchSize = miniBatchSize*10;
end
if ~exist('removeFlatPatches','var')
    removeFlatPatches = false;
end

%% =================== INITIAL E-STEP =================== %%

% extract firstBatchSize patches from the images (columns are patches)
X = dataSource(firstBatchSize);

% if we specify to remove the flat patches than remove any patches where
% the standard deviation of the pixel values is less than 0.002
if removeFlatPatches
    inds = std(X,1)<0.002;
    X(:,inds) = [];
end

% number of sample patches extracted (should be firstBatchSize, right?)
N = size(X,2);

% if the GMM wasn't already initilized then do an initial E-step
if ~isstruct(GMM0)
    
    % randomly sample K patches from the N possible to initlize cluster means
    idx = randsample(N,K);
    m = X(:,idx);
    
    % idenitfy which cluster each patch belongs to
    % TODO: figure out what the subtraction does
    [~,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2),[],1);
    
    % if the number of assigned clusters is different than the number
    % initlized then re-initilize 
    while K ~= length(unique(label))
        idx = randsample(N,K);
        m = X(:,idx);
        [~,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2),[],1);
    end
    
    
    R = full(sparse(1:N,label,1,N,K,N));
    eta = 1;

% otherwise do a nomal E-step
else
    R = zeros(N,K);
    for k = 1:K
        R(:,k) = loggausspdf2(X,GMM.covs(:,:,k))';
    end
    
    R = bsxfun(@plus,R,log(GMM.mixweights));
    T = logsumexp(R,2);
       
    R = bsxfun(@minus,R,T);
    R = exp(R);
    
    eta = (1+T0)^-alpha;
end


%% =================== OPTIMIZE GMM =================== %%

for t=1:nIter
    
    %% =================== M-STEP =================== %%
    s = sum(R,1);
    
    % if there are no zero probabilites there, use this mini batch
%     if (all(s>0))
        GMM.mixweights = GMM.mixweights*(1-eta) + eta*s/N;
        for k = 1:K
            sR = sqrt(R(:,k));
            Xo = bsxfun(@times,X(:,sR>0),sR(sR>0)');
            if s(k)>0
                Sig = double((Xo*Xo')/s(k));
                Sig = Sig + 1e-5*eye(size(Xo,1)); %trace(Sig)/size(Sig,1)*
                [V,D] = eig(Sig);
                D = diag(D);
                D(D<=0) = 1e-5;
                Sig = V*diag(D)*V';
                Sig = (Sig+Sig')/2;
                GMM.covs(:,:,k) = GMM.covs(:,:,k)*(1-eta) + eta*(Sig);
            end
            
        end

    if t<10
        eta = eta/2;
    else
        eta = (t+T0)^-alpha;
    end
    
    %% =================== GET MORE DATA =================== %%
    
    if t<10 %why does he do a different size for the first 10?
        X = dataSource(firstBatchSize);
    else
        X = dataSource(miniBatchSize);
    end
    
    if removeFlatPatches
        inds = std(X,1)<0.002;
        X(:,inds) = [];
    end
    
    %% =================== E-STEP =================== %%
    
    N = size(X,2);
    R = zeros(N,K);

    for k = 1:K

        [V,D] = eigs(GMM.covs(:,:,k),size(GMM.covs,1)-1);
        tt = V'*X;
    
        R(:,k) = -((size(D,1))/2)*log(2*pi) - 0.5*sum(log((diag(D)))) - 0.5*sum(tt.*(D\tt),1)';

    end
    
    R = bsxfun(@plus,R,log(GMM.mixweights));
    T = logsumexp(R,2);
    loglikelihood(t) = sum(T)/N;
    loglikelihood(t) = loglikelihood(t)/(size(X,1)-1)/log(2); % loglikelihood
    
    R = bsxfun(@minus,R,T);
    R = exp(R);
    
    %% =================== PRINT & SAVE INFORMATION =================== %%
    
    if VERBOSE 
        fprintf('Iteration %d of %d, logL: %.2f %s\n',t,nIter,loglikelihood(t),outputFileName);
        subplot(1,2,1);
        plot(loglikelihood(1:t),'o-'); drawnow;
        subplot(1,2,2);
        plot(sort(GMM.mixweights,'descend')); set(gca,'YLim',[0 1]); drawnow;
    end
    
    % save GMM 
    if mod(t,2)==0
        save(outputFileName,'GMM','t','nIter','eta','miniBatchSize','llh','alpha','T0');
    end
    
    
    
end  % end iterations

