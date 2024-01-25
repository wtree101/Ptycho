function [u,mask,errGroup,SNRG,ii]=PhaseRetrievalAP(Y,x,A,At,Params,AtA,init)
%Lasso oversampling without TV 
crop=@(x) x(32:end-32,32:end-32);

flatten=@(x) x(:);


if isfield(Params,'flagBlind')
   flagBlind=Params.flagBlind;
   %disp(["Flag-Blind=",num2str(flagBlind)]);
else 
   flagBlind=0;
end
if flagBlind==1
    cmapidx=Params.cmapidx;
    [Nx,Ny]=size(x);
	Qoverlap=@(frames) reshape(accumarray(flatten(cmapidx),frames(:),[Nx*Ny 1]),Nx,Ny);

end

if ~isfield(Params,'tol')
    tol=eps;
else
    tol=Params.tol;
end

sizeU=size(x);

itmax=Params.itmax;
%Build A matrix
%spectral initialization 

verbose=Params.verbose;

u0=init;

if flagBlind==1
 YMean=sqrt(sum(Y,3)/size(Y,3));
 mask=myifft2((YMean));mask=ifftshift(mask);% approximated mask
%mask = ones(64,64);
A=@(u) myfft2(mask.*u(cmapidx));
maskstack = mask(:,:,ones(size(Y,3),1));
AtA=Qoverlap(abs(maskstack).^2);
end

u=u0;
z0= A(u0);
Ysqrt=sqrt(Y);
%PS=@(z) A(At(z)./AtA);

%PS=@(z) A(real(At(z)./AtA));
%PS=@(z) A((At(z)./AtA));
PM=@(z)(Ysqrt.*sign(z));



reg=max(AtA(:))*1e-3;
if flagBlind==1
    snrC=@(x,x_ref) snrComptBlind(x,x_ref);
else
    snrC=@(x,x_ref) snrComptC(x,x_ref);
end


for ii=1:itmax
       %step I for u
      
       %z=2*beta*PS(PM(z0))+beta*z0-beta*PS(z0)+(1-2*beta)*PM(z0);
	   

	   z=PM(A(u));
      

       

	   
       %% update the variables
       
       %u=(At(z)./AtAS);
	   

	   %u=(At(z)+reg*u)./(AtA+reg);
       u=(At(z))./(AtA);
	   
       %u=(At(z)+reg*u)./AtAS;
       if flagBlind==1
	   	 
		 %BtB=sum(abs(u(cmapidx)),3); %???
         BtB=sum(abs(u(cmapidx)).^2,3); %???square
		 regM=max(BtB(:))*1e-1;
		% mask= (sum(conj(u(cmapidx)).*myifft2(z),3)+regM*mask)./(BtB+regM);
         mask= (sum(conj(u(cmapidx)).*myifft2(z),3))./(BtB);
		 
		 
	     A=@(u) myfft2(mask.*u(cmapidx));
		 At=@(frames) Qoverlap(conj(mask).*myifft2(frames));
		 maskstack = mask(:,:,ones(size(Y,3),1));
         AtA=Qoverlap(abs(maskstack).^2);
		 reg=max(AtA(:))*1e-3;
	   end
	   
	   
       %errGroup(ii)= norm(abs(u(:))-x(:),'fro')/norm(x,'fro');
       errGroup(ii)=sum(vec(abs(abs(A(u))-Ysqrt)))/sum(vec(Ysqrt));
       SNRG(ii)=snrC(crop(u),crop(x));

       if  errGroup(ii)<=tol
        disp(['After',num2str(ii),'_th iteration, AP stopped!'])
        break
       end
       if (verbose==1) & (mod(ii,10)==0) & (flagBlind==0)
       imshow(abs(u),[]);
       title([num2str(ii),'_{th} SNR=',num2str(snrComptC(u,x))]);
       drawnow;
       end
       if (verbose==1) & (mod(ii,10)==0) & (flagBlind==1)
       imshowpair(abs(u),imresize(abs(mask),[size(u,1),size(u,2)]),'montage');
       title([num2str(ii),'_{th} SNR=',num2str(snrComptC(u,x))]);
       drawnow;
       
       end
end


%profile viewer
end


