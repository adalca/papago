function [s, err] = model1sigma(X, W)
% a quick wrapper for sigma model0

    xw = (X .* W);
    sm = xw' * xw;
    s = sm ./ size(X, 1);
% 
%     % method 2. loop.
%     sm = 0;
%     for i = 1:size(X, 1)
%         w = W(i, :)';
%         x = X(i, :)';
%         
%         z = w .* x;
%         q = z * z';
%         sm = sm + q;
%     end
%     s = sm ./ size(X, 1);
     [~, err] = cholcov(s);