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
itmax=2e3;
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
%data_list = ["stacks_random_dist16_blur1.mat","stacks_random_dist8_blur2.mat"];
%data_list = ["stacks_real_random_dist8_sig5_5mat.mat","stacks_real_random_dist8_sig7_0mat.mat","stacks_random_dist8_blur.mat"];
data_list = ["stacks_regular_dist8_sig15_15.mat","stacks_regular_dist8_sig15_15.mat"];
for i=1:size(data_list,2)
    name = data_list(i);
    load(strcat('data\',name)); %read data
    
    
    %MyPathSave=[MyPath,num2str(dist),'-',num2str(jjj)]
    
    a1 = data.image;
    amask = data.phobe;
    u0=ones(size(data.image)); %initial value
    nframes=size(data.stacks,3);
    [Nx,Ny]=size(data.image);
    cmapidx = data.cmapidx;
    %modes = 6;
    modes_list = [12];
    for j=1:size(modes_list,2)
        modes = modes_list(j);
        for flagOrt = 0:1
            %% The simplest AP
            if (modes==1 && flagOrt==1)
                continue;
            end
            
            t0=cputime;
            verbose = 1;
            % ParamsAP
            ParamsAP.tol=TOL;
            ParamsAP.verbose=verbose;
            ParamsAP.flagBlind=flagBlind;
            ParamsAP.flagOrt = flagOrt;
            ParamsAP.cmapidx=cmapidx;
            ParamsAP.itmax=itmax*i; % 2000 and 4000
            ParamsAP.is_project=1;
            %ParamsAP.verbose=1;
            start_time=tic;
            % disp(['Starting to compute PR with AP!'])
            % [uAP,errAP,snrAP,counter,masks] = PhaseRetrievalAP_mixedstates(data.stacks,a1,ParamsAP,ones(size(a1)),1);
            
            output_name = strcat(['outcome/0911AP/APresult_ort',num2str(flagOrt),'modes',num2str(modes),'pr',num2str(i),'iter',num2str(itmax)],name);
            disp(strcat('Starting to compute PR with ADMM!',output_name));
            
            %ParamsAP.InitialPhobe = masks(:,:,1);
            
            
            %uAPs = uAP.* rand(Nx,Ny).*exp(1i*2*pi*rand(Nx,Ny));
            
            
            %ParamsAP.is_project= i;
            tic
            [uAP,Rfactor,snrAP,counter,masks] =PhaseRetrievalAP_mixedstates(data.stacks,a1,ParamsAP,ones(size(a1)),modes);
            % eAP=norm(abs(uAP)-a1,'fro');
            % eAP=norm(abs(uAP-a1),'fro');
            
            result.u = uAP; result.snr = snrAP; result.iter=counter;result.masks=masks;
            %result.e = eAP;
            result.R = Rfactor;
            end_time=toc;
            result.time=end_time-start_time;
            
            
            save(output_name,'result');
            clear result;
            
            
        end
    end
    clear data;
end


return








