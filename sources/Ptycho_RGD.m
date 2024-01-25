function [u,con_ker,errR,snrR,snrR_ker,ii,Rfactor]=Ptycho_RGD(Y,x,Params,cmapidx,amask)
%blind PtychoPR
%global sizeY;global sizeU;

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

[Nx,Ny]=size(x); % x is the complex image 
%cropstack=bsxfun(@times,x(cmapidx),amask);
[px,py,nframes]=size(Y); %px is the size of mask

%% compute \sum_j S^T_j z_j
Qoverlap=@(frames) reshape(accumarray(cmapidx(:),frames(:),[Nx*Ny 1]),Nx,Ny);
flatten=@(x) x(:);
%Qoverlap=@(frames) reshape(accumarray(flatten(cmapidx(:,:,:,ones(modes,1))),frames(:),[Nx*Ny 1]),Nx,Ny);

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


u0=Params.init;

%Initialization for Lambda
if flagGPU==1
    % Lambda=gpuArray(single(zeros(size(z0)))); %GPUArray.zeros(size(z0))
    Lambda=gpuArray.zeros(size(z),'single');
    %    errR=gpuArray(single(zeros(itmax,1)));
    errR=gpuArray.zeros(itmax,1,'single');
    Rfactor=gpuArray.zeros(itmax,1,'single');
else
    %Lambda=zeros(size(z)); 
    errR=(zeros(itmax,1));
    if flagSingle==1
        %Lambda=zeros(size(z),'single');
        %errR=single(zeros(itmax,1));
        errR=zeros(itmax,1,'single');
        Rfactor=zeros(itmax,1,'single');
    end
    
end

snrR=errR+.0;
snrR_ker=errR+.0;


t0=cputime;
if flagGPU==1
    disp('GPU-ADMM starts!')
else
    disp('CPU-ADMM starts!')
    
end

%whos
% Y - sqrt
%u = u0;
u0=Params.init;
Image_recovered = tensor([px*py,Nx*Ny,Nx*Ny],1); % 不知道从0开始好不好。 初始化问题？ 图片都是从全0开始初始化的吗？
ker = fspecial('gaussian', [64 64], 5);
ker_true = fspecial('gaussian', [64 64], 2);
Image_recovered = Image_recovered.initialize(ker,u0);
%Image_recovered_init =cccc Image_recovered;
eta = Params.eta;
measurement = Y.^2;
normalized_mask = reshape(mapminmax(amask(:),0,1),[64 64]);
Qoverlap=@(frames) reshape(accumarray(cmapidx(:),frames(:),[Nx*Nx 1]),Nx,Nx);
for ii=1:itmax
    con_ker = reshape(Image_recovered.U{1},[px,py]);
    v = reshape(Image_recovered.U{2},[Nx,Ny]);
    w = reshape(Image_recovered.U{3},[Nx,Ny]);
    FV = Image_recovered.forward(v,normalized_mask,cmapidx,0);
    FW = Image_recovered.forward(w,normalized_mask,cmapidx,1);
    FVFW = FV.*FW;
    y0 = Image_recovered.A * imfilter(FVFW,con_ker,'circular','conv');

   % y02 = forward(Image_recovered.A,con_ker,v,w,amask,cmapidx);
      

    Grad_F = Image_recovered.get_proj_grad(cmapidx, FV,FW,FVFW,y0 ,measurement,normalized_mask,Qoverlap);
    %disp(sum(Grad_F.U{1}));
    Image_recovered = Image_recovered.retraction(Grad_F, eta);
    u = reshape(Image_recovered.U{2},[Nx,Nx]);

 
    %% Output
    %E=abs(A(u))-Y;
    errR(ii)= norm(u(:)-x(:),'fro')/norm(x,'fro');
    
    %whos  lenY  errR  RF
    
    E=sqrt(abs(y0))-Y; %y0 is square here
    Rfactor(ii)=sum(abs(E(:)))/sum(Y(:));
    
%     snrM=snrComptBlind(mask,amask); snrU=snrComptBlind(u,x);
%     % snrM=0;snrU=0;
%     snrR(ii)=snrU;
    snrU=snrComptBlind(u,x);
    snrR(ii)=snrU;
    snrR_ker(ii) = snrComptBlind(con_ker,ker_true);
    if  errR(ii)<=Params.TOL
        disp(['After',num2str(ii),'_th iteration, ADMM  stopped!'])
        break;
    end
    %
    if (mod(ii,10)==0)
        disp(num2str(ii))
    end
    disp(['Rfactor=',num2str(Rfactor(ii))])
    disp(['snr=',num2str(snrR(ii))])
     disp(['snr_ker=',num2str(snrR_ker(ii))])

    if Params.verbose==1 && (mod(ii,10)==0)
        % colormap(gray)
        
        figure(1);
        %subplot(2,modes,1)
        imshow(abs(u),[]);
        title([num2str(ii),'_{th} SNR=',num2str(snrR(ii))]);
        figure(2);
        %subplot(2,modes,1)
        imshow(abs(reshape(Image_recovered.U{1},[px,px])),[]);

        title([num2str(ii),'_{th} SNR=',num2str(snrR_ker(ii))]);
        figure(3);
        imshow(abs(reshape(Image_recovered.U{3},[Nx,Nx])),[]);
        
        
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



function phi = forward(cof,con_ker,v,w,amask,cmapidx) %target w is conj(v). amask is conj(omega)
    v_stack = v(cmapidx);
    v_stack = v_stack.* amask; %QU？ 
    phiv = myfft2(v_stack);
    w_stack = w(cmapidx);
    w_stack = conj(w_stack) .* amask; %QU？ 
    %phiw = myifft2(w_stack); % forward model has i
    phiw = conj(myfft2(w_stack)); % conj here we use conj on w_stack
    phi = cof * imfilter(phiv .* phiw,con_ker,'circular','conv');

    % should have conjugate.. though I don't think it will make a
    % difference
end



