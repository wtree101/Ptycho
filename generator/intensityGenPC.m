function [Y,kappa,masks]=intensityGenPC(u,amask,cmapidx,Sig)
%myfft2=@(x) fft2(x);

[px,py,~]=size(cmapidx);
A=@(u,m)(myfft2(bsxfun(@times,u(cmapidx), m) ));
%x1=1:nx;x2=1:ny;[X1,X2]=meshgrid(x1,x2);X = [X1(:) X2(:)];


%hsize=30;

%hsize=min(max(Sig(1),Sig(2))*4+1, nx/2-1);
if max(Sig(:))==0
  Y=abs(A(u,amask)).^2;
  disp('Coherent data!')
  kappa=[];
  return
end
hsize=min(Sig*4+1, px/2-1);


ker = fspecial('gaussian',round(hsize),max(Sig));
%ker = fspecial('disk',Sig(1));
%ker = fspecial('average',20);
%kappa=zeros(size(amask));kappa(end/2-hsize/2:end/2+hsize/2-1,end/2-hsize/2:end/2+hsize/2-1)=ker;
kappa=ker;
hsize = size(ker);
%{
kappa= mvnpdf(X, [0,0], [Sig(1) Sig(3); Sig(3) Sig(2)]); % gaussian kernel
kappa = reshape(kappa,length(X1),length(X2));
%}
Y=zeros(size(cmapidx));

t0=tic;
hh=fix(hsize/2)-1;

modes = 4*hh(1)*hh(2);
masks = zeros(px,px,modes);
counter=0;
%for ii=-round(nx/2):1:round(nx/2)-1;    for jj=-round(ny/2):1:round(ny/2)-1;
for ii=-hh(1):hh(1)  
    for jj=-hh(2):hh(2)
        mask=circshift(amask,[ii,jj]);
        factor = ker(ii+hh(1)+1,jj+hh(2)+1);
        Y=Y+abs(A(u,mask)).^2*factor;   
        counter = counter + 1;
        masks(:,:,counter) = mask * sqrt(factor);
    end
end
%{
t0=tic;
for ii=-hsize/2:hsize/2-1;    
    for jj=-hsize/2:hsize/2-1;
        us=circshift(u,[ii,jj]);
        Y=Y+abs(A(us,amask)).^2*kappa(ii+nx/2+1,jj+ny/2+1);       
    end
end
%}

%e=Y1-Y;norm(e(:));Y=Y1;
toc(t0)
YClean=abs(A(u,amask)).^2;
disp(['SNR of data = ', num2str(snrCompt(sqrt(Y(:)),sqrt(YClean(:))))])
end