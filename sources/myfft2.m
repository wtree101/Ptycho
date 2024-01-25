function y = myfft2(x)
%
% normalized fft2
%
[n1,n2,n3] = size(x);
y = fft2(x)/sqrt(n1*n2);
