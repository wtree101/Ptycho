%demo for fully-coherent ptychography  on 7/25/2021

%%
%Pre setting

clear all;
getd = @(p)path(p,path);
MyPath='data/';


%%
% parameters for problem.
% and parameters for solvers
flagBlind=1; 
%noiseFlag=0;
itmax=5e2;
flagSimu=1;
getd('sources/')
if ~exist(MyPath,'dir')
    mkdir(MyPath);
end
verbose=1;
maskType=1;
%imageFlag=2;
TOL=1e-12;

if flagBlind==1
    snrC=@(x,x_ref) snrComptBlind(x,x_ref);
else
    snrC=@(x,x_ref) snrComptC(x,x_ref);
end

%%
%
%load('data/stacks_blur.mat') %read data
load('data/stacks_2_random_dist16.mat')


%MyPathSave=[MyPath,num2str(dist),'-',num2str(jjj)]

a1 = data.image;
amask = data.phobe;
u0=ones(size(data.image)); %initial value
nframes=size(data.stacks,3);
[Nx,Ny]=size(data.image);
cmapidx = data.cmapidx;

%% compute \sum_j S^T_j z_j
Qoverlap=@(frames) reshape(accumarray(cmapidx(:),frames(:),[Nx*Ny 1]),Nx,Ny);
maskstack = data.phobe(:,:,ones(nframes,1));
A=@(u)( myfft2(u(cmapidx).*maskstack)); %A generate stack from origin
At=@(z)(Qoverlap( conj(maskstack).*myifft2(z)) ); % psi to numerator
AtA=Qoverlap(abs(maskstack).^2);  %denominator    
%Y=abs(A(a1)).^2; %isequal(Y,dfstack2) = 1


%% The simplest AP 
% disp(['Starting to compute PR with AP!'])
% t0=cputime;
% 
% % ParamsAP
% ParamsAP.tol=TOL;
% ParamsAP.itmax=itmax;
% ParamsAP.verbose=verbose;
% ParamsAP.flagBlind=flagBlind;
% ParamsAP.cmapidx=cmapidx;
% tic
% [uAP,errAP,snrAP,counter] =PhaseRetrievalAP(data.stacks,a1,  A, At,ParamsAP,AtA,ones(size(a1)));
% 
% toc
% eAP=norm(abs(uAP)-a1,'fro');


%ADMM
disp(['Starting to compute PR with ADMM!'])
t0=cputime;

%
ParamsADMM.itmax=1e2;ParamsADMM.TOL=TOL;
ParamsADMM.verbose=verbose;
ParamsADMM.flagBlind=flagBlind;
ParamsADMM.init=u0;
ParamsADMM.beta=5e-2;
tic
[uADMM,maskADMM,AtA,errADMM,snrADMM]=solverBlindADMMGPU_CXI(data.stacks,a1,ParamsADMM,cmapidx,amask);

toc



return








