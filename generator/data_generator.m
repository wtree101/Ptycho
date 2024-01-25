function data_generator(imageFlag,gridFlag,noiseFlag,blurFlag)
%%
%Pre setting
% data_generator(2,1,0)
%noiseFlag 0 no noise
%gridFlag 1 rectangular
%imageFlag 2 camera man
% blurFlag 1 guassian
%data_generator(2,1,0,1)
%clear all;
getd = @(p)path(p,path);
MyPath='data/';
sizeImg=128;sizeProb=64;%FWHM=16
%flagBlind=1;

ParamsADMM.beta=5e-2;
itmax=1e2;


flagSimu=1;
getd('sources/')
if ~exist(MyPath,'dir')
    mkdir(MyPath);
end

verbose=1;
maskType=1;
%imageFlag=2;
distList=[8];
TOL=1e-12;

% if flagBlind==1
%     snrC=@(x,x_ref) snrComptBlind(x,x_ref);
% else
%     snrC=@(x,x_ref) snrComptC(x,x_ref);
% end

%gridFlag = 1;
dist = 8; %translation

%ini amask
%load('generator\masks\unregular_mask')
[a1,dfstack2,Q,cmapidx,amask,trans]=AuBallsetup_revised(dist,dist,gridFlag,imageFlag,sizeImg,sizeProb,0);

%%
% noise (unfinished)

%
YNoise = dfstack2;
%blur for mixed states (unfinished)
if blurFlag==1 %pc model1
    %YNoise = imgaussfilt(dfstack2,1);
    H = fspecial('gaussian',[64 64],2);
    YNoise = imfilter(dfstack2,H,'circular','conv');
    %YNoise = imgaussfilt(dfstack2,1,'circular','conv');
    %H = fspecial('gaussian',64,1);
    %H = fspecial('average',3);
    %YNoise = imfilter(dfstack2,H,'replicate');
elseif blurFlag==2 %pc model2
    [YNoise,kappa,masks] = intensityGenPC(a1,amask,cmapidx,[15 15]);
    data.masks = masks;
end
%%
if noiseFlag==1
    eta =0.125;
    YNoise=poissrnd(YNoise*(eta))/eta;
    disp(['Noisy intensity SNR=',num2str(snrCompt((YNoise),(dfstack2)))])
end
%% 


data.sizeImg = sizeImg;
data.sizeProb = sizeProb;
data.image = a1;
data.xdist = dist;
data.ydist = dist;
data.nframes = size(dfstack2,3);
%data.ker = H;  
data.cmapidx = cmapidx;
data.phobe = amask;
data.stacks = YNoise;

%%
% save
save('data2/stacks_regular_dist8_clean.mat','data');

end
