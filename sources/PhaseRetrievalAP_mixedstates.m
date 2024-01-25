function [u,errGroup,SNRG,ii,inimask]=PhaseRetrievalAP_mixedstates(Y,x,Params,init,modes)
%Lasso oversampling without TV
% x old image u new one
% 8/6 errGroup - R factor; delete proj on u
% Ort 20/per
crop=@(x) x(32:end-32,32:end-32);

flatten=@(x) x(:);
Ysqrt=sqrt(Y);

cmapidx=Params.cmapidx;
[Nx,Ny]=size(x);
%Qoverlap=@(frames) reshape(accumarray(flatten(cmapidx),frames(:),[Nx*Ny 1]),Nx,Ny);


% multimodes frames: Nx*Ny*nframes*modes

if ~isfield(Params,'tol')
    tol=eps;
else
    tol=Params.tol;
end

if ~isfield(Params,'flagOrt')
    flagOrt=0;
else
    flagOrt=Params.flagOrt;
end

sizeU=size(x);

itmax=Params.itmax;
errGroup=zeros(itmax,1,'single');
errGroup=zeros(itmax,1,'single');
%Build A matrix
%spectral initialization

verbose=Params.verbose;
[px,py,nframes] = size(Y);
u0=init;

% Initialization for mask and u ...
if isfield(Params,'InitialPhobe')
    mask=Params.InitialPhobe;
else
    YMean=sum(Ysqrt,3)/size(Y,3); %why not sqrt
    mask=myifft2((YMean)); mask=ifftshift(mask);% approximated mask
end

%disturb

masks = mask(:,:,ones(modes,1));
maskstacks = mask(:,:,ones(nframes,1),ones(modes,1));
famask = generate_circle(64);
if isfield(Params,'is_project')
    pr = Params.is_project;
else
    pr=1;
end

%amp = mean(mask(:));
rng('default');
for r = 1:modes % randomly genererate
    disturbed = mask .* rand(px,py).*exp(1i*2*pi*rand(px,py));
    masks(:,:,r) = disturbed;
    maskstacks(:,:,:,r) = disturbed(:,:,ones(nframes,1));
end


%Functions
Qoverlap=@(frames) reshape(accumarray(flatten(cmapidx(:,:,:,ones(modes,1))),frames(:),[Nx*Ny 1]),Nx,Ny);
%At=@(frames) Qoverlap(conj(maskstacks).*frames); %denometer frames need
%update?!
%64*64*nframes*modes


u=u0;
%z0= A(u0(cmapidx)); % z - psi
Ysqrt=sqrt(Y);

% if flagBlind==1
%     snrC=@(x,x_ref) snrComptBlind(x,x_ref);
% else
%     snrC=@(x,x_ref) snrComptC(x,x_ref);
% end

snrC=@(x,x_ref) snrComptBlind(x,x_ref);

for ii=1:itmax
    %step I for u
    %z=PM(A(u));
    
    z = Fourier_magnitude_update(modes,Ysqrt,u(cmapidx),maskstacks);
    
    %% update the variables
    
    %u=(At(z)./AtAS);
    
    % object update
    At=@(frames) Qoverlap(conj(maskstacks).*frames); % careful. At needs update!
    AtA=Qoverlap(abs(maskstacks).^2);
    reg=max(AtA(:))*1e-3;
    u=(At(z)+reg*u)./(AtA+reg);
    %u=(At(z))./(AtA); stable terms
    %     index = abs(u)>1;
    %     u(index) = u(index) ./ abs(u(index));
    
    
    %phobes update
    %BtB=sum(abs(u(cmapidx)).^2,3);
    BtB=sum(abs(u(cmapidx)).^2,3); %square?
    regM=max(BtB(:))*1e-3;
    %masks maskstacks 64*64*modes —— mask_stack 64*64*frames*modes
    masks = (squeeze(sum(conj(u(cmapidx)).*z,3)) + regM*masks)./(BtB+regM); % squeeze to drop 1 column
    if (pr==1)
        masks = project(masks,famask);
    end
    %masks = (squeeze(sum(conj(u(cmapidx)).*z,3)))./(BtB);
    %masks= (sum(conj(u(cmapidx)).*z,3)+regM*masks)./(BtB+regM);
    maskstacks = permute(masks(:,:,:,ones(nframes,1)),[1 2 4 3]); %generate masktacks through permute
    %     for k=1:modes
    %         mask = masks(:,:,k);
    %         maskstacks(:,:,:,k) = mask(:,:,ones(nframes,1)); %generate new maskstacks
    %     end
    
    if (mod(ii,20)==0 && flagOrt==1)
        masks = orthogonal(masks);
        maskstacks = permute(masks(:,:,:,ones(nframes,1)),[1 2 4 3]);
    end
    
    
    
    %errGroup(ii)= norm(u(:)-x(:),'fro')/norm(x,'fro'); % ps abs(u) - x wrong as x is complex...
    %errGroup(ii)=sum(vec(abs(abs(A_mixed(u))-Ysqrt)))/sum(vec(Ysqrt));
    zz = A_mixed(u(cmapidx),maskstacks,modes);
    E=sqrt(sum(abs(zz).^2,4))-Ysqrt;
    errGroup(ii)=sum(abs(E(:)))/sum(Ysqrt(:));
    SNRG(ii)=snrC(crop(u),crop(x));
    
    if (mod(ii,10)==0)
        disp(num2str(ii))
    end
    if  errGroup(ii)<=tol
        disp(['After',num2str(ii),'_th iteration, AP stopped!'])
        break
    end
    flagBlind = 1;
    if (verbose==1) & (mod(ii,10)==0) & (flagBlind==0)
        imshow(abs(u));
        title([num2str(ii),'_{th} SNR=',num2str(snrComptC(u,x))]);
        drawnow;
    end
    if (verbose==1) & (mod(ii,10)==0) & (flagBlind==1)
        
        
        % colormap(gray)
        figure(1);
        %subplot(2,modes,1)
        imshow(abs(u),[]);
        title([num2str(ii),'_{th} SNR=',num2str(snrComptC(u,x))]);
        
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
    
    drawnow;
    
end
inimask = masks;
end


%profile viewer


function psi = Fourier_magnitude_update(modes,Ysqrt,u_stack,maskstacks)
%z - psi (64,64,frames,k) maskstacks(64,64,frames,k) u_stack 64*64*frames
%object(64,64,frames,1->k)
% phi 64*64*frames*k
% Ysqrt(64,64,frames)


phi = A_mixed(u_stack,maskstacks,modes); %phi exit wave

psi = phi; %copy for update. don't worry for copy-on-write
amp = sqrt(sum( abs(phi).^2,4)); % 64*64*frames amplitute for updation
for k=1:modes
    
    psi(:,:,:,k) = myifft2( Ysqrt./ amp .* phi(:,:,:,k) );
    %psi(:,:,:,k) = myifft2( new_amp .* phi(:,:,:,k) ); %numerics problems?
    %psi(:,:,:,k) = myifft2( Ysqrt .* sign(phi(:,:,:,k)) );
end

end

function phi = A_mixed(u_stack,maskstacks,modes) %u_stack = u(cmapidx 64*64*nframes)
%generate phi F(OjP)
% u_stack u(cmapidx) u image complex
object = u_stack(:,:,:,ones(modes,1)); %add one dimension
exit_wave = object.*maskstacks;
phi = myfft2(exit_wave);
end

function new_masks = project(masks,famask)
[px,py,modes] = size(masks);
projectors = famask(:,:,ones(1,modes));
fmasks = myfft2(masks);
fmasks = fmasks .* projectors;
new_masks = myifft2(fmasks);
end

% function masks = orthogonal(masks)
% [px,py,modes] = size(masks);
% ss = zeros(px*py,modes);
% for k=1:modes
%     m = masks(:,:,k);
%     ss(:,k) = m(:);
% end
% [q,r]  = qr(ss');
% for k=1:modes
%    masks(:,:,k)=reshape(r(k,:),px,py)
% end
% end

% function masks = orthogonal(masks)
% [px,py,modes] = size(masks);
% ss = zeros(px*py,modes);
% for k=1:modes
%     m = masks(:,:,k);
%     ss(:,k) = m(:);
% end
% [u,s,~] = svd(ss*ss')
% for k=1:modes
%    masks(:,:,k)=reshape(q(:,k),px,py)
% end
% end
% function object_stacks = overlap_projection_objectupdate(psi,maskstacks)
%
% end
%
% function maskstacks = overlap_projection_maskupdate(psi,u)
%
% end


