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
itmax=1e3;
flagSimu=1;
getd('sources/')
if ~exist(MyPath,'dir')
    mkdir(MyPath);
end
verbose=0;
maskType=1;
%imageFlag=2;
TOL=1e-12;
itmax = 3e2;
if flagBlind==1
    snrC=@(x,x_ref) snrComptBlind(x,x_ref);
else
    snrC=@(x,x_ref) snrComptC(x,x_ref);
end

%%
%
%data_list = ["stacks_real_random_dist8_sig5_5mat.mat","stacks_real_random_dist8_sig7_0mat.mat"];
%data_list = ["stacks_random_dist8_blurdisk5.mat","stacks_random_dist8_blurdisk5.mat","stacks_random_dist8_blurdisk1.mat"];
%data_list = ["stacks_real_random_dist8_sig5_5matresult_ADMM300iter_modes9_ort2_2ratio0.0001stacks_real_random_dist8_sig5_5mat.mat.mat","stacks_random_dist8_blur.mat","stacks_real_random_dist8_sig7_0mat.mat"];

%data_list = ["stacks_regular_dist8_blur1_new.mat"]; %,"stacks_random_dist8_blurdisk5.mat"];
%data_list = ["stacks_real_random_dist8_sig15_15_noise_0.0675.mat"];
%data_list = ["stacks_regular_dist8_gu1.mat"];
%"stacks_real_random_dist8_sig15_15_noise_0.25.mat","stacks_real_random_dist8_sig15_15_noise_0.125.mat",...,
%  "stacks_real_random_dist8_sig15_15_noise_0.0675.mat"
% data_list = ["stacks_real_random_dist8_sig15_15_noise_0.25.mat","stacks_real_random_dist8_sig15_15_noise_0.3375.mat",...,
%     "stacks_real_random_dist8_sig15_15_noise_0.01.mat","stacks_real_random_dist8_sig15_15_noise_0.005.mat"];
%data_list = ["stacks_real_random_dist8_sig15_15_noise_0.0675.mat"];
%data_list = ["stacks_real_random_dist8_sig15_15_noise_0.25.mat","stacks_real_random_dist8_sig15_15_noise_0.125.mat","stacks_real_random_dist8_sig15_15_noise_0.03375.mat"];
data_list = ["stacks_regular_dist4_blur2_circ_new.mat"];
%data_list = ["stacks_real_regular_dist8_sig15_15"];
mode_lists{1} = [12];
ratio_lists{1} = [500];
rand_seed = {"default"};
%beta_list = [5 1 5e-1 1e-1 5e-2 1e-2 5e-3 1e-3];
%beta_list = [1e-1 5e-2 1e-2 5e-3 1e-3];
beta_list = [5e-2];
iter_list = ones(1,size(data_list,2));
iter_list(:) = 300;
for i=1:size(data_list,2)
    mode_lists{i} = mode_lists{1};
    ratio_lists{i} = ratio_lists{1};
    beta_lists{i} = beta_list;
end



for i=1:size(data_list,2)
    name = data_list(i);
    itmax = iter_list(i);
    load(strcat('data2/',name)); %read data
    
    %MyPathSave=[MyPath,num2str(dist),'-',num2str(jjj)]
    %load("outcome\result_ADMM_800iter_modes15_ort2_0stacks_random_dist8_blurdisk5.mat");
    
    %% Initialization
    a1 = data.image;
    amask = data.phobe;
    u0=ones(size(data.image)); %initial value
    %u0=result.u;
    nframes=size(data.stacks,3);
    [Nx,Ny]=size(data.image);
    cmapidx = data.cmapidx;
    %modes = 6;
   % flagBlind = 1;
    verbose = 1;
    
    %ini
    %ort_mask = orthogonal(data.masks);
    %ParamsADMM.initmasks = ort_mask(:,:,1:3);
    
    flagBlind = 1;
%     
    %% The ADMM
    %ADMM
    for r_s = 1:size(rand_seed,2)
        t0=cputime;
        %
        ParamsADMM.rand = rand_seed{r_s};
        ParamsADMM.itmax=itmax; %in paper
        ParamsADMM.TOL=TOL;
        ParamsADMM.verbose=verbose;
        ParamsADMM.flagBlind=flagBlind;
        ParamsADMM.init=u0;
       % ParamsADMM.initmasks = data.masks;
        % ParamsADMM.beta=5e-2; %how to set beta?
        ParamsADMM.mode_keep = 6; %select 6 modes
        % ParamsADMM.beta_ort=ParamsADMM.beta*0.2;
        %start_time=tic;
        %  modes_list = [1 3 6 9 15]; % 1 3 6 10 15
        %modes_list = [1 3 6 10 15];
        
        ratio = ratio_lists{i};
        modes_list =  mode_lists{i};
        beta_list =  beta_lists{i};
        %tic
        for b=1:size(beta_list,2)
            for j=1:size(modes_list,2)
                %for flagOrt=0:2:2
                for flagOrt=0:2:0
                    ParamsADMM.ort_iter = 20;
                    if flagOrt==1 || flagOrt==3
                        ParamsADMM.ort_iter=1;
                    end
                    
                    %             if (i==1 && flagOrt==0)
                    %                 continue
                    %             end
                    for r_id = 1:size(ratio,2)
                        if flagOrt~=2 && r_id>1
                            continue
                        end
                        ParamsADMM.beta = beta_list(b);
                        ParamsADMM.beta_ort=ParamsADMM.beta*ratio(r_id);
                        ParamsADMM.flagOrt=flagOrt;
                        modes = modes_list(j);
                        ParamsADMM.initmasks = amask(:,:,ones(modes,1));
                        if (modes==1 && flagOrt==1)
                            continue;
                        end
                        if (flagOrt==3 && modes<= ParamsADMM.mode_keep)
                            continue;
                        end
                        output_name = strcat(['outcome/0524/result_ADMM',num2str(ParamsADMM.itmax),'iter_modes',num2str(modes),'_ort2_',num2str(flagOrt),'ratio', ...,
                            num2str(ratio(r_id)),'beta',num2str(beta_list(b))],'rand',num2str(ParamsADMM.rand),name);
                        disp('Starting to compute PR with ADMM!');
                        disp(output_name);
                        
                        tic
                        [uADMM,Rfactor,snrADMM,counter,masks,cor]=ADMM_mixed_ort_revised(data.stacks,a1,ParamsADMM,cmapidx,amask,modes);
                        %ADMM_mixed(Y,x,Params,cmapidx,amask,modes)
                        end_time=toc;
                        
                        %eADMM=norm(abs(uADMM-a1),'fro');
                        
                        result.u = uADMM;  result.snr = snrADMM; result.iter=counter;result.masks=masks;
                        %result.cor = cor;
                        %result.e = eADMM;
                        result.time = end_time;
                        result.R = Rfactor;
                        %result.time=end_time-start_time;
                        
                        
                        save(output_name,'result');
                        clear result;
                    end
                    
                end
            end
        end
    end
    clear data;
end


return








