function sxr = sx(X, S, dodebug)
% perform S * X, or S(:,:,i) * X(i, :)' for each of i = 1:N
% S is D x D x N
% X is N X D
% sxr is N x D
% trick here: transform S to DxDN matrix, and X to DN x 1;

    [N, D] = size(X);
    assert(all(size(S) == [D, D, N]), 'The size of sigma is incorrect');
    
    % "fast" method
    sxr = multiprod(S, reshape(X', [D, 1, N]));
    sxr = squeeze(sxr)';
    
    % slow/other method
    if nargin == 3 && dodebug
        scheck = zeros(N, D);
        for i = 1:N
            scheck(i, :) = S(:,:,i) * X(i, :)';
        end
        assert(all(isclose(sxr(:), scheck(:))));
    end
end
