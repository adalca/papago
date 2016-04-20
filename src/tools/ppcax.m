function [coeff, score, latent, mu, v, rsltStruct] = ppcax(Y, k, varargin)
%PPCA Probabilistic Principle Component Analysis
%   COEFF = PPCA(Y,K) returns the coefficients for the N by P data matrix Y
%   based on a probabilistic principal component analysis (PPCA) model.
%   Rows of Y correspond to observations and columns to variables. Each
%   column of COEFF contains coefficients for one principal component. The
%   columns are in descending order in terms of component variance. K is
%   the number of principal components requested and it must be a positive
%   scalar integer value and less than the maximum possible rank, min(N,P).
%
%   The PPCA model seeks to relate a P-dimensional observation vector y to
%   a corresponding K-dimensional vector of latent (or unobserved) variable
%   x which is normal with zero mean and covariance of eye(K). The
%   relationship is:
%                           y' = W*x' + mu + err,
%   where the isotropic error term err is Gaussian with zero mean and
%   covariance of v*eye(k) and v is the residual variance. Therefore, under
%   this model, the observations y has mean of mu and covariance of
%   W*W'+v*eye(k). The estimate of W is determined by an Expectation -
%   Maximization (EM) algorithm. When data matrix Y has missing values, the
%   EM algorithm treats the missing values as additional latent variables.
%   At convergence, the columns of W will span the principal subspace, but
%   they are not orthonormal. PPCA obtains COEFF by orthogonalization of W.
%
%   [COEFF,SCORE] = PPCA(Y,K) also returns the principal component score,
%   which is the representation of Y in the principal component space. Rows
%   of SCORE correspond to observations, columns to components.
%
%   [COEFF,SCORE,LATENT] = PPCA(Y,K) also returns the principal component
%   variances.
%
%   [COEFF,SCORE,LATENT,MU] = PPCA(Y,K) also returns the estimated mean of
%   each variable in Y.
%
%   [COEFF,SCORE,LATENT,MU,V] = PPCA(Y,K) returns the isotropic residual
%   variance V.
%
%   [COEFF,SCORE,LATENT,MU,V,S] = PPCA(Y,K) stores final results
%   at convergence in the rsltStruct:
%
%         S.W        - W at convergence.
%         S.Xexp     - Conditional Expectation of X at convergence.
%         S.Recon    - Reconstructed observations using K principle 
%                      components. This is a low-dimension approximation to 
%                      the input data Y, and it is equal to 
%                                 MU + SCORE*COEFF'.
%         S.v        - Residual variance.
%         S.RMSResid - Root mean square of residuals.
%         S.NumIter  - Total number of iteration at convergence.
%         S.nloglk   - Negative log-likelihood at convergence.
%
%   [...] = PPCA(Y,K,'PARAM1',val1,'PARAM2',val2,...) specifies one or more
%   of the following parameter name/value pairs:
%
%          'W0'      - A P-by-K matrix of initial value of W. The default
%                      is a matrix of random values.
%
%          'v0'      - A positive scalar that specifies the initial value 
%                      of v. The default is a random number.
%          
%          'Options' - An options structure as created by the STATSET
%                      function. PPCA uses the following fields:
%             'Display' - Level of display output.  Choices are 'off' (the
%                         default), 'final', and 'iter'.
%             'MaxIter' - Maximum number of steps allowed. The default is
%                         1000. Unlike in optimization settings, reaching
%                         MaxIter is regarded as convergence.
%              'TolFun' - Positive number giving the termination tolerance
%                         for the negative log-likelihood function.  The
%                         default is 1e-6.
%                'TolX' - Positive number giving the convergence threshold
%                         for relative change in the elements of W.
%                         The default is 1e-6.
%
%   Example:
%     >> load hald;
%     >> [coeff,score,latent] = ppca(ingredients,2);
%
%   See also PCA, PCACOV, PCARES, BIPLOT, BARTTEST, CANONCORR, FACTORAN,
%   ROTATEFACTORS.

% Reference: 
%   [1] Tipping, M. E. and Bishop, C. M., Probabilistic Principal Component
%       Analysis. Journal of the Royal Statistical Society: Series B
%       (Statistical Methodology), 1999.
%   [2] Sam Roweis. EM algorithms for PCA and SPCA. In Proceedings of the
%       1997 conference on Advances in neural information processing
%       systems 10 (NIPS '97), 1998, 626-632. 
%   [3] Alexander Ilin and Tapani Raiko. 2010. Practical Approaches to
%       Principal Component Analysis in the Presence of Missing Values. J.
%       Mach. Learn. Res. 11 (August 2010), 1957-2000.

%   Copyright 2012-2014 The MathWorks, Inc.

narginchk(2,inf);

% Validate required input argument Y
if ~(isnumeric(Y)&&min(size(Y))>1)
    error(message('stats:ppca:BadY'));
end

Y = Y';% Y is in PCA format (row observation vector) EM PPCA uses column observation vector

wasNaN = isnan(Y);
hasMissingValue = any(wasNaN(:));
allNaN = all(wasNaN,1);
Y(:,allNaN) = [];
wasNaN(:,allNaN)= [];
obs = ~wasNaN;
numObs = sum(obs(:)); 
[p,n] = size(Y);

% Handling special cases:
if isempty(Y)
    coeff = NaN;
    score = NaN;
    latent = NaN;
    mu = NaN;
    v = NaN;
    rsltStruct = NaN;
    return;
end

% Validate required input argument k
maxRank = min(n,p); % max possible rank of Y
flagWarnK = false;
if ~internal.stats.isScalarInt(k,1)
    error(message('stats:ppca:Badk'));
elseif k > maxRank-1
    k = max(1,maxRank -1);
    flagWarnK = true;
    warning(message('stats:ppca:Warnk',maxRank,k));
end

paramNames = {'Options', 'W0',      'v0'};
dfts =       {[],        randn(p,k), rand};
[Opt,W,v,setFlag] = internal.stats.parseArgs(paramNames,dfts,varargin{:});

% Check 'Options' argument values
defaultopt = statset('ppca');
TolX = statget(Opt,'TolX',defaultopt,'fast');
if ~(isscalar(TolX) && (TolX>0))
    error(message('stats:ppca:BadTol','TolX'));
end
TolFun = statget(Opt,'TolFun',defaultopt,'fast');
if ~(isscalar(TolFun) && (TolFun>0))
    error(message('stats:ppca:BadTol','TolFun'));
end
MaxIter = statget(Opt,'MaxIter',defaultopt,'fast');
DispOpt = statget(Opt,'Display',defaultopt,'fast');
DispOpt = internal.stats.getParamVal(DispOpt,{'off','iter','final'},'Display');
    
% Validate user supplied initial values
if setFlag.W0 && flagWarnK && size(W,2)== maxRank
    W(:,maxRank)=[]; % remove the last column as k is set to maxRank-1
    warning(message('stats:ppca:TruncateW0',n,k));
end
if setFlag.W0 && ...
        (~(isnumeric(W)&& ismatrix(W) && isequal(size(W),[p,k])) || any(isnan(W(:))))
    error(message('stats:ppca:BadW0',p,k));
end
if setFlag.v0 && (~(isscalar(v)&&v>0) || isnan(v) || v==inf)
    error(message('stats:ppca:Badv0'));
end

% Suppress undesired warnings.
warnState = warning('query','all');
warning('off','MATLAB:illConditionedMatrix');
warning('off','MATLAB:nearlySingularMatrix');
cleanupObj = onCleanup(@() warning(warnState));

% Preallocate memory
mu = zeros(p,1);
X = zeros(k,n);
Wnew = zeros(p,k);
C = zeros(k,k,n);
nloglk = inf;

dispnum = strcmp(DispOpt,{'off','iter','final'});

headernames = {'Iteration'    'Variance'   '|Delta X|' 'Negative Log-likelihood'};
header = sprintf('\n%9s %12s %12s   %12s:\n', headernames{:});
iterfmtstr = '%9d %12f %12f   %12f\n';
if dispnum(2)
    fprintf(header);
end

itercount = 0;

if hasMissingValue
    % If Y has any missing value, use the following algorithm
    while (itercount < MaxIter)
        itercount = itercount +1;
        
        for j = 1:n
            y = Y(:,j);
            idxObs = obs(:,j);
            w = W(idxObs,:);
            
            % Use Sherman-Morrison formula to find the inv(v.*eye(k)+w'*w)
            Cj = eye(k)/v-w'*w/(eye(k)+w'*w/v)/(v^2);
            X(:,j) = Cj*(w'*(y(idxObs)-mu(idxObs)));
            C(:,:,j) = Cj;
        end
        
        % update mean
        mu = nanmean(Y-W*X,2);
        
        % update W
        for i = 1:p
            idxObs = obs(i,:);
            M = X(:,idxObs)*X(:,idxObs)'+v.*sum(C(:,:,idxObs),3);
            ww = X(:,idxObs)*(Y(i,idxObs)-mu(i))';
            Wnew(i,:) = M\ww;
        end
        
        % update residual variance
        vsum = 0;
        for j=1:n
            wnew = Wnew(obs(:,j),:);
            vsum = vsum + sum( ( ...
                (Y(obs(:,j),j) - wnew*X(:,j) - mu(obs(:,j))).^2 + v*diag(wnew*C(:,:,j)*wnew') ) );
        end
        vnew = vsum/numObs;
        
        % Compute negative log-likelihood function
        % nloglk_new = ppcax_incomplete_nlogl(obs, Y, mu, Wnew, vnew);
        nloglk_new = -inf;
        
        dw = max(max(abs(W-Wnew) / (sqrt(eps)+max(max(abs(Wnew))))));
        
        W = Wnew;
        v = vnew;
        
        if dw<TolX
            break;
        elseif (nloglk-nloglk_new) < TolFun
            break;
        end

        nloglk = nloglk_new;

        if dispnum(2)
            fprintf(iterfmtstr,itercount, v,dw,nloglk);
        end
        
        if itercount == MaxIter
            warning(message('stats:ppca:MaxIterReached', MaxIter));
        end
    end % End of While Loop
    
    % Center X
    muX = mean(X,2);
    X = bsxfun(@minus,X,muX);
    % Update the mean of Y, mu
    mu = mu + W*muX;
    
else % if no missing value, use the faster algorithm for complete data
    [W,X,mu,v,itercount,dw,nloglk]...
        = emppca_complete(Y,k,W,v,MaxIter,TolFun,TolX,dispnum,iterfmtstr);
    if all(W(:)==0)
        coeff = zeros(p,k);
        coeff(1:(p+1):end)=1;
        score = zeros(n,k);
        latent = zeros(k,1);
        mu = mean(Y,2)';
        v = 0;        
        rsltStruct.W = W;
        rsltStruct.Xexp = X'; % make it consistent with the input data format
        rsltStruct.Recon = repmat(mu,n,1); % add estimated mean back
        rsltStruct.v = v;
        rsltStruct.NumIter = itercount;
        rsltStruct.RMSResid = 0;
        rsltStruct.nloglk = nloglk;
        return;
    end
end
% Reconstruction:
WTW = W'*W;
Y_hat = W/(WTW)*(WTW+v*eye(k))*X; % Centered, will add bias term back later

% Compute the RMS Residuals
diff = bsxfun(@minus,Y-Y_hat,mu);
diff(~obs)=0;
rmsResid = norm(diff,'fro')/sqrt(numObs);

% Print out the final results at convergence
if ~dispnum(1)
    if ~dispnum(2)
        fprintf(header);
        fprintf(iterfmtstr,itercount, v,dw,nloglk);
    end
    fprintf('\n%s.\n',getString(message('stats:ppca:FinalRMSResidual',...
        sprintf('%g',rmsResid))));
end

% Orthogonalize W to the standard PCA subspace
[coeff,~] = svd(W,'econ');
score = Y_hat'*coeff;
[~,latent]=eig(score'*score);
latent = sort(diag(latent),'descend')/(n-1);
mu = mu';

% Enforce a sign convention on the coefficients -- the largest element in
% each column will have a positive sign.
[~,maxind] = max(abs(coeff), [], 1);
[d1, d2] = size(coeff);
colsign = sign(coeff(maxind + (0:d1:(d2-1)*d1)));
coeff = bsxfun(@times, coeff, colsign);
if nargout > 1
    score = bsxfun(@times, score, colsign); % scores = score
end

% Insert NaNs back to rows that has all NaN elements
score = internal.stats.insertnan(allNaN,score);

% Store additional information at convergence
if nargout > 5
    rsltStruct.W = W;
    rsltStruct.Xexp = X'; % make it consistent with the input data format
    rsltStruct.Recon = bsxfun(@plus,Y_hat',mu); % add estimated mean back
    rsltStruct.v = v;
    rsltStruct.NumIter = itercount;
    rsltStruct.RMSResid = rmsResid;
    rsltStruct.nloglk = nloglk;
end

smallLatent = latent < (max(size(Y))*eps(max(latent)));
if any(smallLatent)
    error(message('stats:ppca:Singular',sum(smallLatent),...
        sum(~smallLatent)));
end
% --------------End of Main Function---------------------------------------


% ------------- Local Subfunction------------------------------------------
function [Wnew,Xmu,mu,vnew,iter,dw,nloglk_new] = emppca_complete(Y,k,W,v,MaxIter,TolFun,TolX,dispnum,fmtstr)
% local function. EM algorithm for the complete data case

[p,n] = size(Y); % p - dimension; n - number of observations
mu = mean(Y,2); % mean is the sample mean
Y = bsxfun(@minus,Y,mu); % Centering the data

iter = 0;
nloglk = inf;

traceS = Y(:)'*Y(:)/(n-1);

while iter < MaxIter
    iter = iter +1;
    SW = Y*(Y'*W)/(n-1);
    M = W'*W+v*eye(k);
    
    Wnew = SW/(v*eye(k)+M\W'*SW);
    vnew = (traceS-trace(SW/M*Wnew'))/p;

    % Check Convergence:
    dw = max(max(abs(W-Wnew) / (sqrt(eps)+max(max(abs(Wnew))))));
    dv = abs(v-vnew)/(eps+v);
    delta = max(dw,dv);
    
    CC = Wnew*Wnew'+vnew*eye(p);
    nloglk_new = (p*log(2*pi) + log(det(CC)) + trace(CC\Y*Y'/(n-1)) )*n/2;
    
    W = Wnew;
    v = vnew;
    
    if dispnum(2)
        fprintf(fmtstr,iter,v,delta,nloglk_new);
    end
    
    if delta < TolX
        break;
    elseif (nloglk - nloglk_new) < TolFun
        break;
    elseif abs(vnew) < sqrt(eps)
        break;
    end

    nloglk = nloglk_new;
      
    if iter == MaxIter
        warning(message('stats:ppca:MaxIterReached', MaxIter));
    end
end

Xmu = M\Wnew'*Y;