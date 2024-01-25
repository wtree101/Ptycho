function [xx,yy,rr]=make_grid(img,ctr,shift,gpu)
if ~exist('gpu')
  gpu=0;
end


  
[M,N]=size(img);

if M*N==1
    M=img;N=img;
elseif M*N==2
    M=img(2);N=img(1);
end
% if ~gpu
%   [xx,yy]=meshgrid(1:N,1:M);
% else
%   [xx,yy]=meshgrid(gsingle(1:N),gsingle(1:M));
% end
%[xx,yy]=meshgrid(gsingle(1:N),gsingle(1:M));
[xx,yy]=meshgrid(1:N,1:M);
  
if nargin>1 %centering
    if ctr
        xx=xx-mean(xx(1,:)+1/2);
        yy=yy-mean(yy(:,2)+1/2);
    end
%     if ctr==.5
% %    xx=floor(xx);
% ^    yy=floor(yy);    
end
if nargin==3
    if shift
    xx=fftshift(xx);
    yy=fftshift(yy);
    end
end

    
if nargout>2 %rr dist from center
    rr=sqrt(xx.^2+yy.^2);
end

  
