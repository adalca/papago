function idx = lenMatrix2idx(lenMat)
% go from a vector of lengths to a cell of indices
%
% e.g. [4, 8, 2] --> {1:4, 5:12, 13:14}

    cumsumR = cumsum(lenMat(:)');

    rends = cumsumR; 
    rstarts = [1 cumsumR(1:end-1) + 1];

    idx = reshape(arrayfunc(@(a,b) a:b, rstarts, rends), size(lenMat));