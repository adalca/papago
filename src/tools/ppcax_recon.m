function Y_hat = ppcax_recon(Y, W, v, mu)

    Y = Y';
    mu = mu';
    obs = ~isnan(Y);
    k = size(W, 2);

    
    X = zeros(k, size(Y,2));
    for j = 1:size(Y,2)
        y = Y(:,j);
        idxObs = obs(:,j);
        w = W(idxObs,:);

        % Use Sherman-Morrison formula to find the inv(v.*eye(k)+w'*w)
        Cj = eye(k)/v-w'*w/(eye(k)+w'*w/v)/(v^2);
        X(:,j) = Cj*(w'*(y(idxObs)-mu(idxObs)));
    end
    
    WTW = W'*W;
    Y_hat = W/(WTW)*(WTW+v*eye(k))*X; % Centered, will add bias term back later
    
    Y_hat = bsxfun(@plus, Y_hat', mu');
    
