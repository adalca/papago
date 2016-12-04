function wgs = wv2sigma(wgs)

    for k = 1:size(wgs.params.W, 3)
        wgs.params.sigma(:,:,k) = wgs.params.W(:,:,k) * wgs.params.W(:,:,k)' + ...
            wgs.params.sigmasq(k) * eye(size(wgs.params.W, 1));
    end
