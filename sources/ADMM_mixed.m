function [u,errR,snrR,ii,masks]=ADMM_mixed(Y,x,Params,cmapidx,amask,modes)
%blind PtychoPR
%global sizeY;global sizeU;

cpu2GPU=@(x) gpuArray(single(x));
if ~isfield(Params,'flagGPU')
    flagGPU=0;flagSingle=0;
else
    flagGPU=Params.flagGPU;
    flagSingle=1;
end

if ~isfield(Params,'flagOrt')
    flagOrt=0;
else
    flagOrt=Params.flagOrt;
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

% if Params.verbose==1
%     figure;
% end

TOL=Params.TOL;

[NxM,NyM]=size(amask);
crop=@(x)x(fix(NxM/2):end-fix(NxM/2),fix(NyM/2):end-fix(NyM/2));

[Nx,Ny]=size(x);
%cropstack=bsxfun(@times,x(cmapidx),amask);
[px,py,nframes]=size(Y);

%% compute \sum_j S^T_j z_j
%Qoverlap=@(frames) reshape(accumarray(cmapidx(:),frames(:),[Nx*Ny 1]),Nx,Ny);
flatten=@(x) x(:);
Qoverlap=@(frames) reshape(accumarray(flatten(cmapidx(:,:,:,ones(modes,1))),frames(:),[Nx*Ny 1]),Nx,Ny);

itmax=Params.itmax;
beta=Params.beta;

%if flagGPU==1

%Y=single(int16(Y));
%% Initializing

%Initializtion for mask
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

masks = mask(:,:,ones(modes,1));
maskstacks = mask(:,:,ones(nframes,1),ones(modes,1));

%amp = mean(mask(:));
for r = 1:modes % randomly genererate
    disturbed = mask .* rand(px,py).*exp(1i*2*pi*rand(px,py));
    masks(:,:,r) = disturbed;
    maskstacks(:,:,:,r) = disturbed(:,:,ones(nframes,1));
end

%A=@(u)( myfft2(u(cmapidx).*maskstacks));
%At=@(z)(Qoverlap( conj(maskstacks).*myifft2(z)) );
%A=@(u)( myfft2(u(cmapidx).*maskstack));At=@(z)(Qoverlap( conj(maskstack).*myifft2(z)) );
%AtA=Qoverlap(abs(maskstacks).^2);
%figure(2);imshow(abs(mask),[]);colorbar;
u0=Params.init;
z0=A_mixed(u0(cmapidx),maskstacks,modes);% let |z0|\leq a %A_mixed(u_stack,maskstacks,modes)

%Initialization for Lambda
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
        masks0=masks;
        maskAV=sum(abs(u0(cmapidx)).^2,3);
        epsM=max(maskAV(:))*1e-5;
        %zeroRegM=(abs(maskAV)<=eps);    maskAV=maskAV.*(1-zeroRegM)+zeroRegM;
        %zhatFinv=myifft2(z0+Lambda/beta);
        %zhatFinv=myifft2(z0);
        %mask=(sum(conj(u0(cmapidx)).*zhatFinv,3))./maskAV.*(1-zeroRegM)+mask0.*zeroRegM;
        masks=( (squeeze(sum(conj(u0(cmapidx)).*myifft2(z0),3))) + epsM*masks0 )./(maskAV+epsM);
        maskstacks = permute(masks(:,:,:,ones(nframes,1)),[1 2 4 3]); %generate masktacks through permute
    end
    if (mod(ii,20)==0)
        masks = orthogonal(masks);
        maskstacks = permute(masks(:,:,:,ones(nframes,1)),[1 2 4 3]);
    end
    
    %% step II
    
    % A=@(u)( myfft2(u(cmapidx).*maskstacks));
    At=@(zFinv)(Qoverlap( zFinv.*conj(maskstacks)) );%reduce duplicative cost
    AtA=Qoverlap(abs(maskstacks).^2);
    
    epsU=max(AtA(:))*1e-3;
    %zeroRegU=(abs(AtA)<=eps);
    %AtA=AtA.*(1-zeroRegU)+zeroRegU;
    u=(At(myifft2(z0))+epsU*u0)./(AtA+epsU); %u new u-2
    
    %% step III
    
    zz=A_mixed(u(cmapidx),maskstacks,modes);
    %tmp=A(u)-Lambda/beta;
    z0=zz-Lambda; % A - lambda
    %tmpflag=abs(tmp)~=0;
    
    %    z=(Y+beta*abs(tmp))/(1+beta).*(sign(tmp).*tmpflag)+(1-tmpflag).*z0;
    %z0=(Y+beta*abs(z0))/(1+beta).*sign(z0);
    amp = sqrt(sum( abs(z0).^2,4)); % 64*64*frames amplitute for updation
    for k=1:modes
        z0(:,:,:,k) =  ( (Y ./ amp).*z0(:,:,:,k) +  beta*z0(:,:,:,k) ) / (1+beta); %???
    end
    
    %% update the variables
    
    time(ii)=cputime-t0;
    
    
    u0=u;
    %Lambda=Lambda+beta*ee;
    % updata multipliers
    Lambda= Lambda + (z0-zz); %zo-zz z0 new z zz A?
    
    %% Output
    %E=abs(A(u))-Y;
    errR(ii)= norm(abs(u(:))-x(:),'fro')/norm(x,'fro');
    %whos  lenY  errR  RF
    %errR(ii)=sum(abs(E(:)))/sum(Y(:));
    
    snrM=snrComptBlind(mask,amask); snrU=snrComptBlind(u,x);
    % snrM=0;snrU=0;
    snrR(ii)=snrU;
    
    if  errR(ii)<=Params.TOL
        disp(['After',num2str(ii),'_th iteration, ADMM  stopped!'])
        break;
    end
    %
    
    if Params.verbose==1 && (mod(ii,10)==0)
        % colormap(gray)
        figure(1);
        %subplot(2,modes,1)
        imshow(abs(u),[]);
        title([num2str(ii),'_{th} SNR=',num2str(snrR(ii))]);
        
        figure(2)
        %colormap(gray)
        for k=1:modes
            subplot(1,modes,k)
            imshow(abs(masks(:,:,k)),[]);
            title([num2str(k),'th mode'])
            axis equal
            axis tight
        end
        
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

function phi = A_mixed(u_stack,maskstacks,modes) %u_stack = u(cmapidx 64*64*nframes)
%generate phi F(OjP)
% u_stack u(cmapidx) u image complex
object = u_stack(:,:,:,ones(modes,1)); %add one dimension
exit_wave = object.*maskstacks;
phi = myfft2(exit_wave);
end

