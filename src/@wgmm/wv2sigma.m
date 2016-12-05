function sigma = wv2sigma(wgs)

    [dHigh, ~, K] = size(wgs.params.W);

    sigma = zeros(dHigh, dHigh, K);
    for k = 1:size(wgs.params.W, 3)
        sigma(:,:,k) = wgs.params.W(:,:,k) * wgs.params.W(:,:,k)' + ...
            wgs.params.sigmasq(k) * eye(size(wgs.params.W, 1));
    end
