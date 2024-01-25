function B = cropmat(A,bsize)
%
% usage: B = cropmat(A,bsize)
%

[n,m] = size(A);
if (n < bsize(1) || m < bsize(2))
   fprintf('image size must be larger than %d\n', bsize);
   exit(1); 
end

imin = floor((n-bsize(1))/2)+1;
jmin = floor((m-bsize(2))/2)+1;
B = A(imin:imin+bsize(1)-1,jmin:jmin+bsize(2)-1);
