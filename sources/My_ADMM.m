function [u,mask,AtA,errR,snrR]=solverBlindADMMGPU_CXI(Y,x,Params,cmapidx,amask)
%blind PtychoPR
%global sizeY;global sizeU;

cpu2GPU=@(x) gpuArray(single(x));
if ~isfield(Params,'flagGPU')
flagGPU=0;flagSingle=0;
else
flagGPU=Params.flagGPU;
flagSingle=1;
end




if isfield(Params,'flagBlind')
   flagBlind=Params.flagBlind;
   %disp(["Flag-Blind=",num2str(flagBlind)]);
else 
   flagBlind=0;
end
if flagBlind==0
mask=amask;
else

YMean=sqrt(mean(Y,3));
mask=myifft2((YMean));mask=ifftshift(mask);% approximated mask

end

if flagGPU==1
    Y=sqrt(cpu2GPU(Y));
    x=cpu2GPU(x);
    %cmapidx=int16(cmapidx);
    amask=cpu2GPU(amask);
    cmapidx=cpu2GPU(cmapidx);
else
    Y=sqrt(Y);
end

TOL=Params.TOL;

[NxM,NyM]=size(amask);
crop=@(x)x(fix(NxM/2):end-fix(NxM/2),fix(NyM/2):end-fix(NyM/2));

[Nx,Ny]=size(x);
%cropstack=bsxfun(@times,x(cmapidx),amask);
nframes=size(Y,3);

%% compute \sum_j S^T_j z_j
Qoverlap=@(frames) reshape(accumarray(cmapidx(:),frames(:),[Nx*Ny 1]),Nx,Ny);

itmax=Params.itmax;
beta=Params.beta;


%if flagGPU==1

%Y=single(int16(Y));


%% testing codes

%%





maskstack = mask(:,:,ones(nframes,1));
A=@(u)( myfft2(u(cmapidx).*maskstack));At=@(z)(Qoverlap( conj(maskstack).*myifft2(z)) );
%A=@(u)( myfft2(u(cmapidx).*maskstack));At=@(z)(Qoverlap( conj(maskstack).*myifft2(z)) );


AtA=Qoverlap(abs(maskstack).^2);

%figure(2);imshow(abs(mask),[]);colorbar;



u0=Params.init;


z0=sign(A(u0)).*(Y);% let |z0|\leq a

%% initialization by ADMM

%u0=real(At(z0)./AtA);
%u0=abs(At(z0)./AtA);
%u0=(At(z0)./AtA);
%z0= A(u0);




if Params.verbose==1
    figure;
end


if flagGPU==1
    % Lambda=gpuArray(single(zeros(size(z0)))); %GPUArray.zeros(size(z0))
    Lambda=gpuArray.zeros(size(z0),'single');
    %    errR=gpuArray(single(zeros(itmax,1)));
    errR=gpuArray.zeros(itmax,1,'single');
else
    Lambda=zeros(size(z0)); errR=(zeros(itmax,1));
    if flagSingle==1
        Lambda=zeros(size(z0),'single');
        %errR=single(zeros(itmax,1));
        errR=zeros(itmax,1,'single');
    end
    
end

snrR=errR+.0;

%Lam=Lambda/beta;
% init mask
%
%mask=randn(size(amask))+1i*randn(size(amask));
%u0=x;z0=A(u0);



%[beta beta2 tau ]
%flagFI
t0=cputime;
if flagGPU==1
    disp('GPU-ADMM starts!')
else
    disp('CPU-ADMM starts!')
    
end

%whos
% Y - sqrt
for ii=1:itmax
    
    z0=z0+Lambda;
    
    if flagBlind==1
    %% step I for illumination with mask0
    mask0=mask;
    maskAV=sum(abs(u0(cmapidx)),3);
    epsM=max(maskAV(:))*1e-5;
    %zeroRegM=(abs(maskAV)<=eps);    maskAV=maskAV.*(1-zeroRegM)+zeroRegM;
    %zhatFinv=myifft2(z0+Lambda/beta);
    
    
    
    %zhatFinv=myifft2(z0);
    %mask=(sum(conj(u0(cmapidx)).*zhatFinv,3))./maskAV.*(1-zeroRegM)+mask0.*zeroRegM;
    mask=(sum(conj(u0(cmapidx)).*myifft2(z0),3)+epsM*mask0)./(maskAV+epsM);
    end
    
    
    
    %% step II
    
    maskstack = mask(:,:,ones(nframes,1));
    A=@(u)( myfft2(u(cmapidx).*maskstack));
    Bt=@(zFinv)(Qoverlap( zFinv.*conj(maskstack)) );%reduce duplicative cost
    AtA=Qoverlap(abs(maskstack).^2);
    
    epsU=max(AtA(:))*1e-3;
    %zeroRegU=(abs(AtA)<=eps);
    %AtA=AtA.*(1-zeroRegU)+zeroRegU;
    u=(Bt(myifft2(z0))+epsU*u0)./(AtA+epsU);
    
    
    zz=A(u);
    
    %tmp=A(u)-Lambda/beta;
    z0=zz-Lambda; % A - lambda
    
    
    
    %tmpflag=abs(tmp)~=0;
    
    %    z=(Y+beta*abs(tmp))/(1+beta).*(sign(tmp).*tmpflag)+(1-tmpflag).*z0;
    z0=(Y+beta*abs(z0))/(1+beta).*sign(z0);
    
    %% update the variables
    
    time(ii)=cputime-t0;
    
    u0=u;
    
    %Lambda=Lambda+beta*ee;
    Lambda=Lambda+ (z0-zz); %zo-zz z0 new z zz A?
    
    
    
    
    E=abs(A(u))-Y;
    %whos  lenY  errR  RF
    errR(ii)=sum(abs(E(:)))/sum(Y(:));

    snrM=snrComptBlind(mask,amask); snrU=snrComptBlind(u,x);
    % snrM=0;snrU=0;
    snrR(ii)=snrU;
	
	
	if  errR(ii)<=Params.TOL
        disp(['After',num2str(ii),'_th iteration, ADMM  stopped!'])
        break;
    end
    %
	
    if Params.verbose==1 && (mod(ii,10)==0)
        %whos
        
        title(num2str(ii))
        subplot(121);
        imshow(abs(crop(u)),[]);
        title([num2str(ii),'_{th} SNR=',num2str(snrU)]);
        drawnow;
        subplot(122);
        imshow(abs(mask),[]);title([num2str(snrM)]);
        drawnow;
        
        disp([num2str(ii),'_{th} iteration--ADMM'])
        
    end
    
    
    
    
    
    
    
    
end


errR=errR(1:ii);
snrR=snrR(1:ii);

if flagGPU==1
u=gather(u);
mask=gather(mask);
AtA=gather(AtA);
errR=gather(errR);
snrR=gather(snrR);
end
%whos


end

