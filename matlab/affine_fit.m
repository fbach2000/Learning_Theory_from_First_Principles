function [a,b] = affine_fit(x,y);

% making sure we have column vectors
if size(x,1)==1, 
    x = x';
end
if size(y,1)==1, 
    y = y';
end
n = length(x);
x = [ x, ones(n,1) ];
 
M = (x'*x) \ (x'*y);
a = M(1);
b = M(2);

