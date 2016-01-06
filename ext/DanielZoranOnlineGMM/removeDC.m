function [X,ns,xs,ys] = removeDC(X,ns,xs,ys)

X = X - repmat(mean(X,1), [size(X,1) 1]); 
