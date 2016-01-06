function s = iwAiw(W, A)
% for each i=1:N performs.
%   s = WW' \hadamard A   (i.e. (W * W') .* A)
%   which is the same as: s = D_W A D_W (D_W is diag(W))
% for W of size NxD and A of size DxD. 
%
% This gets complicated since DxDxN can get memory intensive, which WW' would be in wgmm

    % for size(W, 1) == 1, this is simple:
    % s = W * W' .* A; % fast but requires huge memory
    % s = diag(W) * A * diag(W); % slowest, requires huge memory due to diag.
    % s = bsxfun(@times, bsxfun(@times, W, A), W'); % middle-speed, low memory overhead.
    
    [~, D] = size(W);
    assert(all(size(A, 1) == D));
    assert(all(size(A, 2) == D));
    
    % we'll use the third implementation
    b1 = bsxfun(@times, permute(W, [2, 3, 1]), A);
    s = bsxfun(@times, b1, permute(W, [3, 2, 1]));
end
