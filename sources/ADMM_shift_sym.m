function [u, errR,snrR,ii,masks,cor]=ADMM_shift_sym(Y,x,Params,cmapidx,amask,modes)
%blind PtychoPR
%global sizeY;global sizeU;
% % 8/6 errGroup - R factor; delete proj on u

cpu2GPU=@(x) gpuArray(single(x));
if ~isfield(Params,'flagGPU')
    flagGPU=0;flagSingle=0;
else
    flagGPU=Params.flagGPU;
    flagSingle=1;
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

if ~isfield(Params,'flagOrt')
    flagOrt=0;
else
    flagOrt=Params.flagOrt;
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
beta_ort = Params.beta_ort;
b_ratio = beta_ort/beta;
%if flagGPU==1

%Y=single(int16(Y));
%% Initializing;

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
    %mask=myifft2(ifftshift(YMean));
    % Y 值在四周哦 —— z也是哦
    mask=myifft2(YMean);
    mask=fftshift(mask);% approximated mask
end

%fre domain
famask = generate_circle(64,1);
tamask = generate_circle(64,3.4);
famask = ifftshift(famask);

%Lambda 
lambda = ones(modes,1);
%shifts = zeros(modes,2);

if isfield(Params,'initmasks')
    masks = Params.initmasks;
    maskstacks = permute(masks(:,:,:,ones(nframes,1),ones(2,1)),[1 2 4 3 5]);
    %disp(["Flag-Blind=",num2str(flagBlind)]);
else
% inimask   
    %amp = mean(mask(:));
    % if flagBlind==1
    
    if ~isfield(Params,'rand') 
        rng('default'); %don forget!!
    else
        rng(Params.rand)
    end
    %shifts = floor(randn(modes,2)*(px/2)/4);
    %shifts(1,1) = 0; shifts(1,2) = 0;
    shifts = [0 0;0 30;30 0;30 30;30 -30];
    masks = mask(:,:,ones(modes,1),ones(2,1));
    %maskstacks = mask(:,:,ones(nframes,1),ones(modes,1));
    lambda(1,1) = 1;
    for r = 1:modes % randomly genererateb
        masks(:,:,r,1) = circshift(mask,[shifts(r,1),shifts(r,2)]) * lambda(r);
        masks(:,:,r,2) = circshift(mask,[-shifts(r,1),-shifts(r,2)]) * lambda(r);
    end
    maskstacks = permute(masks(:,:,:,:,ones(nframes,1)),[1 2 5 3 4]);
  
end


%% ini u
u0=Params.init;

u = u0;
ustack = u(cmapidx);
ustacks = ustack(:,:,:,ones(modes,1),ones(2,1));
ustacks_lambda = ustack(:,:,:,ones(modes,1),ones(2,1)); % *lambda
%ini z
zz=A_mixed(ustack,maskstacks,modes);% let |z0|\leq a %A_mixed(u_stack,maskstacks,modes)
z = zz;


%Initialization for ADMM
if flagGPU==1
    % Lambda=gpuArray(single(zeros(size(z0)))); %GPUArray.zeros(size(z0))
    Lambda=gpuArray.zeros(size(z),'single');
    %    errR=gpuArray(single(zeros(itmax,1)));
    
    errR=gpuArray.zeros(itmax,1,'single');
else
    Lambda=zeros(size(z)); errR=(zeros(itmax,1));
    if flagSingle==1
        Lambda=zeros(size(z),'single');
        %errR=single(zeros(itmax,1));
        errR=zeros(itmax,1,'single');
        cor = zeros(itmax,1,'single');
    end
    
end

%Ort
snrR=errR+.0;

options = optimoptions(@fminunc,'MaxIterations',2,'display','off');

%% begin ADMM
t0=cputime;
if flagGPU==1
    disp('GPU-ADMM starts!')
else
    disp('CPU-ADMM starts!')
    
end

%whos
% Y - sqrt



for ii=1:itmax
    %% step III update z
    %tmp=A(u)-Lambda/beta;
    z0=zz-Lambda; % A - lambda
    %tmpflag=abs(tmp)~=0;
    %    z=(Y+beta*abs(tmp))/(1+beta).*(sign(tmp).*tmpflag)+(1-tmpflag).*z0;
    %z0=(Y+beta*abs(z0))/(1+beta).*sign(z0);
    amp = sqrt(sum( abs(z0).^2,4)); % 64*64*frames amplitute for updation
    for k=1:modes
        z(:,:,:,k) =  ( (Y ./ amp).*z0(:,:,:,k) +  beta*z0(:,:,:,k) ) / (1+beta); %???
    end
    
    %new z
    
    z0=z+Lambda;
    Finv_z0 = myifft2(z0);
    %% step II update u
    % A=@(u)( myfft2(u(cmapidx).*maskstacks));
    At=@(zFinv)(Qoverlap( zFinv.*conj(maskstacks(:,:,:,:,1)+maskstacks(:,:,:,:,2))) );%reduce duplicative cost
    AtA=Qoverlap(abs(maskstacks(:,:,:,:,1)).^2 + abs(maskstacks(:,:,:,:,2)).^2);
    epsU=max(AtA(:))*1e-3;
    %zeroRegU=(abs(AtA)<=eps);
    %AtA=AtA.*(1-zeroRegU)+zeroRegU;
    u=(At(Finv_z0)+epsU*u)./(AtA+epsU); %u new u-2
    %     index = abs(u)>1;
    %     u(index) = u(index) ./ abs(u(index)); wrong!
     %% step I for illumination with mask0
    if flagBlind==1
       %% update mask
        ustack = u(cmapidx);
        ustacks_d4 = ustack(:,:,:,ones(modes,1));
        
        for k=1:modes
            ustacks(:,:,:,k,1) = circshift(ustacks_d4(:,:,:,k),[-shifts(k,1) -shifts(k,2)]);
            ustacks(:,:,:,k,2) = circshift(ustacks_d4(:,:,:,k),[shifts(k,1) shifts(k,2)]);
            ustacks_lambda(:,:,:,k,1) = lambda(k) * ustacks(:,:,:,k,1);
            ustacks_lambda(:,:,:,k,2) = lambda(k) * ustacks(:,:,:,k,2);
        end
            
        
        masks0=mask;
        maskAV=squeeze(sum(abs(ustacks_lambda).^2,[3,4,5])); %squeeze?
        epsM=max(maskAV(:))*1e-5;
        %zeroRegM=(abs(maskAV)<=eps);    maskAV=maskAV.*(1-zeroRegM)+zeroRegM;
        %zhatFinv=myifft2(z0+Lambda/beta);
        %zhatFinv=myifft2(z0);
        %mask=(sum(conj(u0(cmapidx)).*zhatFinv,3))./maskAV.*(1-zeroRegM)+mask0.*zeroRegM;
%         if flagOrt==2 % with ort
%             omega0 = omega + Lambda2;
%             masks=( (squeeze(sum(conj(u(cmapidx)).*myifft2(z0),3))) + epsM*masks0 + b_ratio*omega0) ...,
%                 ./(maskAV+epsM+b_ratio);
%         else
        mask=( squeeze(sum( conj(ustacks_lambda).*Finv_z0,[3,4,5])) + epsM*masks0 )./(maskAV+epsM);
        mask = project(mask,famask,tamask);
        
        %% update lambda
        aAV   = squeeze(sum((abs(ustacks).^2) .* (abs(mask).^2),[1,2,3,5]));
       
        %lambda = squeeze(sum( real(conj(Finv_z0).*(ustacks.*mask)),[1,2,3]))./ aAV; % ? squeeze?
        %lambda = max(lambda,0);
        lambda = squeeze(sum( Finv_z0.*conj(ustacks.*mask),[1,2,3,5]))./ aAV;
        abs(lambda)
%         %% update shifts
%         obj = @(shifts) objective(Finv_z0,mask,ustacks_d4,lambda,shifts);
%         [shifts_fr2,~] = fminunc(obj,shifts(2:end,:),options);
%         shifts(2:end,:) = floor(shifts_fr2);
        %% masks update finally
        for r = 1:modes % randomly genererateb
            masks(:,:,r,1) = circshift(mask,[shifts(r,1),shifts(r,2)]) * lambda(r);
            masks(:,:,r,2) = circshift(mask,[-shifts(r,1),-shifts(r,2)]) * lambda(r);
        end
       maskstacks = permute(masks(:,:,:,ones(nframes,1),ones(2,1)),[1 2 4 3 5]); %generate masktacks through permute
       
    end
    
   
    
   
    %% update the multiplier variables
    
    zz=A_mixed(ustack,maskstacks,modes);
    
    Lambda = Lambda + (z-zz); %zo-zz z0 new z zz A?
    
  
    %% Output
    %E=abs(A(u))-Y;
    % errR(ii)= norm(u(:)-x(:),'fro')/norm(x,'fro');
    
    %whos  lenY  errR  RF
    
    E=sqrt(sum(abs(zz).^2,4))-Y;
    errR(ii)=sum(abs(E(:)))/sum(Y(:));
    cor(ii) = coherence(masks);
    
    snrM=snrComptBlind(mask,amask); snrU=snrComptBlind(u,x);
    % snrM=0;snrU=0;
    snrR(ii)=snrU;
    
    if  errR(ii)<=Params.TOL
        disp(['After',num2str(ii),'_th iteration, ADMM  stopped!'])
        break;
    end
    %
    if (mod(ii,10)==0)
        disp(num2str(ii))
    end
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
exit_wave = object.*(maskstacks(:,:,:,:,1) + maskstacks(:,:,:,:,2));
phi = myfft2(exit_wave);
end


function new_mask = project(mask,famask,tamask)

fmask = myfft2(mask);
fmask = fmask .* famask;
new_mask = myifft2(fmask);
new_mask =  new_mask .* tamask;
end

function obj = objective(Finv_z,mask,ustacks_d4,lambda,shifts_fr2)
[~,~,nframes,modes] = size(ustacks_d4);
masks = mask(:,:,ones(modes,1));
masks(:,:,1) = mask * lambda(1);
shifts_fr2 = floor(shifts_fr2);
%maskstacks(:,:,:,1) = mask_r(:,:,ones(nframes,1));
for r = 2:modes % randomly genererate
    mask_r = circshift(mask,[shifts_fr2(r-1,1),shifts_fr2(r-1,2)]) * lambda(r);
    masks(:,:,r) = mask_r;
    %maskstacks(:,:,:,r) = mask_r(:,:,ones(nframes,1));
end
maskstacks = permute(masks(:,:,:,ones(nframes,1)),[1 2 4 3]);
obj = sum(abs(Finv_z - ustacks_d4 .* maskstacks).^2, 'all');
end 







