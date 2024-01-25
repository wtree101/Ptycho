 %demo for fully-coherent ptychography  on 7/25/2021

%%
%Pre setting

clear all;
getd = @(p)path(p,path);
MyPath='data/';


%%
% parameters for problem.
% and parameters for solvers
flagBlind=0;
%noiseFlag=0;
%itmax= 2e4;  %3e3;
itmax = 5e4;
flagSimu=1;
getd('sources/')
if ~exist(MyPath,'dir')
    mkdir(MyPath);
end
verbose=0;
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
data_list = ["stacks_regular_dist8_blur2_circ_new.mat"];
%data_list = ["stacks_regular_dist8_h20.mat"];
%data_list = ["stacks_random_dist8_blurdisk5.mat","stacks_random_dist8_blurdisk5.mat","stacks_random_dist8_blurdisk1.mat"];
%data_list = ["stacks_real_random_dist8_sig5_5mat.mat","stacks_real_random_dist8_sig7_0mat.mat"];
% cores = feature('numcores');
% parpool('local',cores);

for i=1:size(data_list,2)
    name = data_list(i);
    load(strcat('data2/',name)); %read data
    
    %MyPathSave=[MyPath,num2str(dist),'-',num2str(jjj)]
    %load("outcome\result_ADMM_800iter_modes15_ort2_0stacks_random_dist8_blurdisk5.mat");
    
    a1 = data.image;
    amask = data.phobe;
    u0=ones(size(data.image)); %initial value
    %u0=result.u;
    nframes=size(data.stacks,3);
    [Nx,Ny]=size(data.image);
    cmapidx = data.cmapidx;
    %modes = 6;
    
    %% The ADMM
    
    %flagBlind = 1;
    %verbose = 1;
    %ADMM
    
    t0=cputime;
    
    %
    ParamsRGD.itmax=itmax; %in paper
    ParamsRGD.TOL=TOL;
    ParamsRGD.verbose=verbose;
    ParamsRGD.flagBlind=0;
    ParamsRGD.init=u0;
    %ParamsADMM.initmasks = result.masks;
    ParamsRGD.beta=5e-2; %how to set beta?
    ParamsRGD.flagGPU = 0;
    cores = feature('numcores');
    parpool('local',cores);
    %start_time=tic;
  %  modes_list = [1 3 6 9 15]; % 1 3 6 10 15
    %modes_list = [1 3 6 10 15];
    eta_list = [1e-10];
    %tic
    for j=1:size(eta_list,2)
        
        eta = eta_list(j); ParamsRGD.eta = eta;
        disp(['Starting to compute PR with RGD!','eta',num2str(eta)]);
        disp(name);
        %start_time=tic;
        tic;
        %[uADMM,errADMM,snrADMM,counter,masks, Rfactor]=Ptycho_RGD(data.stacks,a1,ParamsRGD,cmapidx,amask,modes);
        %ADMM_mixed(Y,x,Params,cmapidx,amask,modes)
        %data.stacks - Y; data.u - a1; 
        [uRGD,con_ker,errRGD,snrRGD,snr_kerRGD,counter,Rfactor]=Ptycho_RGD_GPU(data.stacks,a1,ParamsRGD,cmapidx,amask);
        %[uRGD,errRGD,snrRGD,counter,Rfactor]=Ptycho_RGD_f(data.stacks,a1,ParamsRGD,cmapidx,amask);
        end_time=toc
        %toc
        
        eRGD=norm(abs(uRGD-a1),'fro');
        
        result.u = uRGD; result.err = errRGD; result.snr = snrRGD; result.iter=counter;%result.masks=masks;
        result.e = eRGD; result.time = end_time;
        result.snr_ker = snr_kerRGD;
        result.R = Rfactor;
        result.ker = con_ker;
        result.time=end_time;
        disp(num2str(result.time))
        output_name = strcat(['outcome/0630/result_RGD_blur2_ss',num2str(ParamsRGD.itmax),'iter'],name);
        save(output_name,'result');
        clear result;
        
        
    end
    clear data;
end


return
exit;








