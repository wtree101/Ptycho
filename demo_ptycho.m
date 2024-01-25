%demo for fully-coherent ptychography  on 7/25/2021

%%
%Pre setting

clear all;
getd = @(p)path(p,path);
MyPath='data/';
sizeImg=128;sizeProb=64;%FWHM=16
flagBlind=1; 
noiseFlag=0; 
ParamsADMM.beta=5e-2;
itmax=1e2;


flagSimu=1;
getd('sources/')
if ~exist(MyPath,'dir')
    mkdir(MyPath);
end
verbose=1;
maskType=1;
imageFlag=2;
distList=[16];
TOL=1e-12;

if flagBlind==1
    snrC=@(x,x_ref) snrComptBlind(x,x_ref);
else
    snrC=@(x,x_ref) snrComptC(x,x_ref);
end

%%
%

for jjj=1
    gridFlag=jjj; %girdflag=1?
    for iii=1:length(distList) % dist = 16 dist?
        [jjj,iii]
		
        dist=distList(iii);
        MyPathSave=[MyPath,num2str(dist),'-',num2str(jjj)]
        
        %generate stacks
        [a1,dfstack2,Q,cmapidx,amask,trans]=AuBallsetup(dist,dist,gridFlag,imageFlag,sizeImg,sizeProb);
        
        u0=ones(size(a1));
        nframes=size(dfstack2,3);
        [Nx,Ny]=size(a1);
        %% compute \sum_j S^T_j z_j
        Qoverlap=@(frames) reshape(accumarray(cmapidx(:),frames(:),[Nx*Ny 1]),Nx,Ny);
        maskstack = amask(:,:,ones(nframes,1));
        A=@(u)( myfft2(u(cmapidx).*maskstack)); %A generate stack from origin
        At=@(z)(Qoverlap( conj(maskstack).*myifft2(z)) ); % psi to numerator
        AtA=Qoverlap(abs(maskstack).^2);  %denominator    
        Y=abs(A(a1)).^2; %isequal(Y,dfstack2) = 1
        
        
        
        %noise setting
        if noiseFlag==1
            YNoise=poissrnd(Y*eta.^2)/eta.^2;
            disp(['Noisy intensity SNR=',num2str(snrCompt((YNoise),(Y)))])
        else
            YNoise=Y;
        end
        
        sizeU=size(a1);
        maskAV=sum(abs(a1(cmapidx)).^2,3);
        [max(AtA(:))/min(AtA(:)),max(maskAV(:))/min(maskAV(:))]
        disp(['Sampling Measurements Ratio is ', num2str(numel(Y)/numel(a1))]);
        
        
        
        Params.flagBlind=flagBlind;
        Params.dist=dist;
        Params.itmax=itmax;
        Params.TOL=TOL;
        Params.noiseFlag=noiseFlag;

        Params.verbose=verbose;
        Params.init=u0;
        Params.flagSimu=flagSimu;
        Params.trans=trans;
		
		

        %% The simplest AP 
        disp(['Starting to compute PR with AP!'])
        t0=cputime;
		ParamsAP.tol=TOL;
		ParamsAP.itmax=itmax;
		ParamsAP.verbose=verbose;
		ParamsAP.flagBlind=flagBlind;
		ParamsAP.cmapidx=cmapidx;
        tic
        [uAP,errAP,snrAP,counter] =PhaseRetrievalAP(YNoise,a1,  A, At,ParamsAP,AtA,ones(size(a1)));
        
        toc
        eAP=norm(abs(uAP)-a1,'fro');
        
        
        %% ADMM
        disp(['Starting to compute PR with ADMM!'])
        t0=cputime;
		
		ParamsADMM.itmax=itmax;ParamsADMM.TOL=TOL;
        ParamsADMM.verbose=verbose;
        ParamsADMM.flagBlind=flagBlind;
		ParamsADMM.init=u0;
        tic
		[uADMM,maskADMM,AtA,errADMM,snrADMM]=solverBlindADMMGPU_CXI(YNoise,a1,ParamsADMM,cmapidx,amask);
		
		toc
		
				

    end
end

return








