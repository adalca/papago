function dst = wtdst(ZI, ZJ, patchSize)

    nVox = prod(patchSize(1:3));

    % extract data and coordinates.
    dataidx = 1:nVox;
    locidx = nVox*2+1:size(ZI, 2);
    X = ZI(:, [dataidx, locidx]);
    Y = ZJ(:, [dataidx, locidx]);

    % weights
    wtidx = nVox + (1:nVox);
    singlewt = ZI(:, wtidx);
    multiwt = ZJ(:, wtidx);
    wt = bsxfun(@times, singlewt, multiwt);
    wt = [wt, ones(size(wt, 1), numel(locidx))];
    wt = bsxfun(@rdivide, wt, sum(wt, 2));


    % apply weights, do sum and take sqrt
    % squares of differences (includes both data and weighted coordinates)
    m = bsxfun(@minus, X, Y);
    ms = m .^ 2;
    mult = wt .* ms;
    s = sum(mult, 2); 
    dst = sqrt(s);