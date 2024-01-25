function rf=rfactor(cmapidx,u,mask,Y,bkg)
 YSqrt=sqrt(Y);
 if ~exist('bkg','var')
 res=abs(myfft2(mask.*u(cmapidx)))-YSqrt;
 else
 res=sqrt(abs(myfft2(mask.*u(cmapidx))).^2+bkg.^2)-YSqrt;
 
 end
 rf=sum(abs(res(:)))/sum(YSqrt(:));
end
