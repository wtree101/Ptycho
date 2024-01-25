function y = myifft2(x)
%
% normalized ifft2
%
[n1,n2,n3]=size(x);
y = ifft2(x)*sqrt(n1*n2);
