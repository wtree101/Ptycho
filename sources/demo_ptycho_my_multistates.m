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
load('data/stacks_blur.mat') %read data


%MyPathSave=[MyPath,num2str(dist),'-',num2str(jjj)]

a1 = data.image;
amask = data.phobe;
u0=ones(size(data.image)); %initial value
nframes=size(data.stacks,3);
[Nx,Ny]=size(data.image);
cmapidx = data.cmapidx;
modes = 5;

%% The simplest AP 
disp(['Starting to compute PR with AP!'])
t0=cputime;

% ParamsAP
ParamsAP.tol=TOL;
ParamsAP.itmax=itmax;
ParamsAP.verbose=verbose;
ParamsAP.flagBlind=flagBlind;
ParamsAP.cmapidx=cmapidx;
tic
[uAP,errAP,snrAP,counter] =PhaseRetrievalAP_mixedstates(data.stacks,a1,ParamsAP,ones(size(a1)),modes);

toc
eAP=norm(abs(uAP)-a1,'fro');






return








