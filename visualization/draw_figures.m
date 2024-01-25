%modes and comparison with AP
name_list = ["AP,r=12","ADMM,r=1","ADMM,r=3","ADMM r=6","ADMM,r=12"]
result_list = ["APresultnew_ort0modes12stacks_random_dist8_blur.mat", ...,
    "result_ADMM_200iter_modes1_ort0stacks_random_dist8_blur.mat", ...,
    "result_ADMM_200iter_modes3_ort0stacks_random_dist8_blur", ...,
    "result_ADMM_200iter_modes6_ort0stacks_random_dist8_blur", ...,
    "result_ADMM_200iter_modes12_ort0stacks_random_dist8_blur"]
mydraw_compare_1(result_list,name_list)
 % Ort experiment
name_list = ["Ort","Not ort"]
result_list = ["result_ADMM_200iter_modes6_ort1stacks_random_dist8_blur.mat", ...,
    "result_ADMM_200iter_modes6_ort0stacks_random_dist8_blur", ...,
   ]

%modes
%density compare

mydraw_compair_2("result_ADMM_200iter_modes12_ort0stacks_random_dist8_blur","model1_gu1_pro1.mat")

mydraw_compair_2("result_ADMM300iter_modes15_ort2_0stacks_random_dist8_model2_sig5_5mat.mat","model1_gu1_pro1.mat"
%ex2
name_list = ["r=1,(7 0)","r=3 (7 0)","r=3 (5 5)","r=15 (5 5)"];
result_list = ["result_ADMM300iter_modes3_ort2_0stacks_random_dist8_model2_sig7_0mat.mat", ...
    "result_ADMM300iter_modes15_ort2_0stacks_random_dist8_model2_sig7_0mat.mat", ...
    "result_ADMM300iter_modes3_ort2_0stacks_random_dist8_model2_sig5_5mat.mat" ,...
    "result_ADMM300iter_modes15_ort2_0stacks_random_dist8_model2_sig5_5mat.mat"];


name_list = ["r=1,(7 0)","r=3 (7 0)","r=3 (5 5)","r=15 (5 5)"];
result_list = ["result_ADMM300iter_modes1_ort2_0stacks_random_dist8_model2_sig7_0mat.mat", ...
    "result_ADMM300iter_modes3_ort2_0stacks_random_dist8_model2_sig7_0mat.mat", ...
    "result_ADMM300iter_modes3_ort2_0stacks_random_dist8_model2_sig5_5mat.mat" ,...
    "result_ADMM300iter_modes15_ort2_0stacks_random_dist8_model2_sig5_5mat.mat"];





mydraw_compare_2("result_ADMM300iter_modes15_ort2_0stacks_random_dist8_model2_sig5_5mat.mat",data.masks
%%%%%%%%%%%
name_list = ["Ort","Not ort"]
result_list = ["APresultnew_ort1modes3stacks_random_dist8_blur.mat","APresultnew_ort0modes3stacks_random_dist8_blur.mat"]
%%
name_list = ["ratio=1","ratio=0.1","ratio=0.01","ratio=0.001","ratio=0.0001","no Ort"];
result_list = ["0808\result_ADMM300iter_modes9_ort2_2ratio1stacks_real_random_dist8_sig5_5mat.mat", ...
    "0808\result_ADMM300iter_modes9_ort2_2ratio0.1stacks_real_random_dist8_sig5_5mat.mat", ...
    "0808\result_ADMM300iter_modes9_ort2_2ratio0.01stacks_real_random_dist8_sig5_5mat.mat" ,...
    "0808\result_ADMM300iter_modes9_ort2_2ratio0.001stacks_real_random_dist8_sig5_5mat.mat",...
     "0808\result_ADMM300iter_modes9_ort2_2ratio0.0001stacks_real_random_dist8_sig5_5mat.mat",...
     "0807\result_ADMM300iter_modes9_ort2_0stacks_real_random_dist8_sig5_5mat.mat"];
 
 
 name_list = ["ratio=0.01","no Ort"];
result_list = [
    "0808\result_ADMM300iter_modes9_ort2_2ratio0.01stacks_real_random_dist8_sig5_5mat.mat" ,...
     "0807\result_ADMM300iter_modes9_ort2_0stacks_real_random_dist8_sig5_5mat.mat"];
%%
name_list = ["ratio=0.001","no Ort"];
result_list = [
    "0808\result_ADMM300iter_modes9_ort2_2ratio0.01stacks_random_dist8_blur.mat" ,...
     "0807\result_ADMM300iter_modes9_ort2_0stacks_random_dist8_blur.mat"];
 
 
 name_list = ["ratio=0.001","ratio=0.0001","no Ort"];
 result_list = [
    "0808\result_ADMM300iter_modes6_ort2_2ratio0.001stacks_real_random_dist8_sig5_5mat.mat" ,...
     "0808\result_ADMM300iter_modes6_ort2_2ratio0.0001stacks_real_random_dist8_sig5_5mat.mat",...
     "0807\result_ADMM300iter_modes6_ort2_0stacks_real_random_dist8_sig5_5mat.mat"];
 
 
 name_list = ["ratio=0.001","ratio=0.0001","no Ort"];
 result_list = [
    "0808\result_ADMM300iter_modes6_ort2_2ratio0.001stacks_random_dist8_blur.mat" ,...
     "0808\result_ADMM300iter_modes6_ort2_2ratio0.0001stacks_random_dist8_blur.mat",...
     "0807\result_ADMM300iter_modes6_ort2_0stacks_random_dist8_blur.mat"];

 %%
 name_list = ["ratio=0.001","ratio=1","No Ort"];
 result_list = [
    "0808\result_ADMM1000iter_modes12_ort2_2ratio0.001stacks_random_dist8_blurdisk5.mat" ,...
     "0808\result_ADMM1000iter_modes12_ort2_2ratio1stacks_random_dist8_blurdisk5.mat",...
     "0809\result_ADMM_800iter_modes15_ort2_0stacks_random_dist8_blurdisk5.mat"];
    
 %%
 [~,int2] = orthogonal_full(st_masks);
plot(cumsum(int2))
hold on 
plot(cumsum(intensities))
load("data\stacks_real_random_dist8_sig10_10mat.mat")
masks_st = data.masks;
[~,int2] = orthogonal_full(masks_st);
 
s = cumsum(int2)
%%
 name_list = ["ratio=0.001","ratio=1","Ort 0","Ort 1"];
 result_list = [
    "0809\result_ADMM300iter_modes6_ort2_2ratio0.001stacks_real_random_dist8_sig10_10mat.mat" ,...
     "0809\result_ADMM300iter_modes6_ort2_2ratio1stacks_real_random_dist8_sig10_10mat.mat",...
     "0809\result_ADMM300iter_modes6_ort2_0ratio0.001stacks_real_random_dist8_sig10_10mat.mat",...
     "0809\result_ADMM300iter_modes6_ort2_1ratio0.001stacks_real_random_dist8_sig10_10mat.mat"];
 
 %% ort ratio=1 ! success. 
 name_list = ["ratio=1","ratio=0.1","Ort 0"];
result_list = [
"0809_2\result_ADMM300iter_modes6_ort2_2ratio1stacks_real_random_dist8_sig10_10mat.mat" ,"0809_2\result_ADMM300iter_modes6_ort2_2ratio0.1stacks_real_random_dist8_sig10_10mat.mat",...,
"0809\result_ADMM300iter_modes6_ort2_0ratio0.001stacks_real_random_dist8_sig10_10mat.mat",]
mydraw_compare_1(result_list,name_list)

%% ratio=1 seems good
name_list = ["ratio=1","Ort 0"];
result_list = [
"0809_2\result_ADMM300iter_modes9_ort2_2ratio1stacks_real_random_dist8_sig10_10mat.mat" ,...,
"0809\result_ADMM300iter_modes9_ort2_0ratio0.001stacks_real_random_dist8_sig10_10mat.mat",]
mydraw_compare_1(result_list,name_list)
%%
name_list = ["ratio=1","ratio=0.6","ratio=0.2"];
result_list = [
"0809_2\result_ADMM300iter_modes12_ort2_2ratio1stacks_real_random_dist8_sig15_15mat.mat" ,...,
"0809_2\result_ADMM300iter_modes12_ort2_2ratio0.6stacks_real_random_dist8_sig15_15mat.mat",...,
"0809_2\result_ADMM300iter_modes12_ort2_2ratio0.2stacks_real_random_dist8_sig15_15mat.mat"];
mydraw_compare_1(result_list,name_list)

%%
% ort
name_list = ["ratio=50","ratio=10","ratio=1","ratio=0.6","ratio=0.2","Not ort"];
result_list = [
"0810\result_ADMM300iter_modes6_ort2_2ratio50stacks_real_random_dist8_sig15_15mat.mat" ,...,
"0810\result_ADMM300iter_modes6_ort2_2ratio10stacks_real_random_dist8_sig15_15mat.mat",...,
"0809_2\result_ADMM300iter_modes6_ort2_2ratio1stacks_real_random_dist8_sig15_15mat.mat",...,
"0809_2\result_ADMM300iter_modes6_ort2_2ratio0.6stacks_real_random_dist8_sig15_15mat.mat",...,
"0809_2\result_ADMM300iter_modes6_ort2_2ratio0.2stacks_real_random_dist8_sig15_15mat.mat",...,
"0810\result_ADMM300iter_modes6_ort2_0ratio50stacks_real_random_dist8_sig15_15mat.mat"];
mydraw_compare_1(result_list,name_list)



%% no effects
name_list = ["ratio=1","ratio=0.6","No ort"];
result_list = [
"0809_2\result_ADMM1200iter_modes15_ort2_2ratio1stacks_random_dist8_blurdisk5.mat" ,"0809_2\result_ADMM1200iter_modes15_ort2_2ratio0.6stacks_random_dist8_blurdisk5.mat",...,
"0809_2\result_ADMM1200iter_modes15_ort2_0ratio1stacks_random_dist8_blurdisk5.mat",]
mydraw_compare_1(result_list,name_list)

%% no effects
name_list = ["ratio=10","ratio=1","No ort"];
result_list = [
"0810\result_ADMM300iter_modes3_ort2_2ratio10stacks_random_dist8_blur.mat" ,"0810\result_ADMM300iter_modes3_ort2_2ratio1stacks_random_dist8_blur.mat",...,
"0810\result_ADMM300iter_modes3_ort2_0ratio1stacks_random_dist8_blur.mat",]
mydraw_compare_1(result_list,name_list)

%% rng default
%% no effectcs for ort  %example,report
% ADMM is faster
name_list = ["ratio=10","ratio=1","No ort","Ort 1","AP Ort0", "AP Ort 1"];
result_list = [
"0810_2\result_ADMM300iter_modes6_ort2_2ratio10stacks_random_dist8_blur.mat" ,...,
"0810_2\result_ADMM300iter_modes6_ort2_2ratio1stacks_random_dist8_blur.mat",...,
"0810_2\result_ADMM300iter_modes6_ort2_0ratio1stacks_random_dist8_blur.mat",...,
"0810_2\result_ADMM300iter_modes6_ort2_1ratio1stacks_random_dist8_blur.mat",...,
"0806\APresult_ort0modes6stacks_random_dist8_blur.mat",...,
"0806\APresult_ort1modes6stacks_random_dist8_blur.mat"
]
mydraw_compare_1(result_list,name_list)
%%
name_list = ["ratio=10","ratio=1","No ort","Ort 1","AP Ort0", "AP Ort 1"];
result_list = [
"0810_2\result_ADMM300iter_modes6_ort2_2ratio10stacks_random_dist8_blur.mat" ,...,
"0810_2\result_ADMM300iter_modes6_ort2_2ratio1stacks_random_dist8_blur.mat",...,
"0810_2\result_ADMM300iter_modes6_ort2_0ratio1stacks_random_dist8_blur.mat",...,
"0810_2\result_ADMM300iter_modes6_ort2_1ratio1stacks_random_dist8_blur.mat",...,
"0806\APresult_ort0modes6stacks_random_dist8_blur.mat",...,
"0806\APresult_ort1modes6stacks_random_dist8_blur.mat"
]
mydraw_compare_1(result_list,name_list)

%% ort effects report
%15 15 and blur1
% some effects modes=9
name_list = ["ratio=1000","raio=100","raio=10","no ort","ort 1","Ort1 1"];
 result_list = [
    "0812\result_ADMM300iter_modes6_ort2_2ratio1000stacks_real_random_dist8_sig15_15mat.mat" ,...
     "0812\result_ADMM300iter_modes6_ort2_2ratio100stacks_real_random_dist8_sig15_15mat.mat",...
     "0812\result_ADMM300iter_modes6_ort2_2ratio10stacks_real_random_dist8_sig15_15mat.mat",...
     "0812\result_ADMM300iter_modes6_ort2_0ratio100stacks_real_random_dist8_sig15_15mat.mat",...
     "0812\result_ADMM300iter_modes6_ort2_1ratio100stacks_real_random_dist8_sig15_15mat.mat"];
mydraw_compare_1(result_list,name_list)

%modes=6 not so good
name_list = ["ratio=1000","raio=100","raio=10","no ort","ort 1"];
 result_list = [
    "0812\result_ADMM300iter_modes6_ort2_2ratio1000stacks_real_random_dist8_sig15_15mat.mat" ,...
     "0812\result_ADMM300iter_modes6_ort2_2ratio100stacks_real_random_dist8_sig15_15mat.mat",...
     "0812\result_ADMM300iter_modes6_ort2_2ratio10stacks_real_random_dist8_sig15_15mat.mat",...
     "0812\result_ADMM300iter_modes6_ort2_0ratio100stacks_real_random_dist8_sig15_15mat.mat",...
     "0812\result_ADMM300iter_modes6_ort2_1ratio100stacks_real_random_dist8_sig15_15mat.mat"];
mydraw_compare_1(result_list,name_list)
%% blur1 
name_list = ["raio=500","raio=250","ratio=100","ratio=10","no ort","ort 1"];
 result_list = [
   
     "0812\result_ADMM300iter_modes6_ort2_2ratio500stacks_random_dist8_blur.mat",...
     "0812\result_ADMM300iter_modes6_ort2_2ratio250stacks_random_dist8_blur.mat",...
     "0812\result_ADMM300iter_modes6_ort2_2ratio100stacks_random_dist8_blur.mat",...         
     "0812\result_ADMM300iter_modes6_ort2_2ratio10stacks_random_dist8_blur.mat",...
      "0812\result_ADMM300iter_modes6_ort2_0ratio100stacks_random_dist8_blur.mat",...
       "0812\result_ADMM300iter_modes6_ort2_1ratio100stacks_random_dist8_blur.mat"];
mydraw_compare_1(result_list,name_list)
%% modes experiment new %report clean
name_list = ["r=6","r=9","r=12","AP,r=12"];
 result_list = [
   
     "0812\result_ADMM300iter_modes6_ort2_0ratio100stacks_random_dist8_blur.mat",...
     "0812\result_ADMM300iter_modes9_ort2_0ratio100stacks_random_dist8_blur.mat",...
     "0812\result_ADMM300iter_modes12_ort2_0ratio100stacks_random_dist8_blur.mat",...
     "0806\APresult_ort0modes12stacks_random_dist8_blur.mat"
    ];
mydraw_compare_1(result_list,name_list)

%%
%noise 2 re
name_list = ["r=12,0","r=12,3","r=9,0","r=9,3"];
 result_list = [
    "0810_2\result_ADMM300iter_modes12_ort2_0ratio1stacks_real_random_dist8_blur_noise_0.0025.mat" ,...
     "0810_2\result_ADMM300iter_modes12_ort2_3ratio1stacks_real_random_dist8_blur_noise_0.0025.mat",...
     "0810_2\result_ADMM300iter_modes9_ort2_0ratio1stacks_real_random_dist8_blur_noise_0.0025.mat",...
     "0810_2\result_ADMM300iter_modes9_ort2_3ratio1stacks_real_random_dist8_blur_noise_0.0025.mat"];
mydraw_compare_1(result_list,name_list)

%some effects but all fails 0.0025 re
name_list = ["r=12,0","r=12,3","r=9,0","r=9,3"];
result_list = [
"0810_2\result_ADMM300iter_modes12_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.0025.mat" ,...
"0810_2\result_ADMM300iter_modes12_ort2_3ratio1stacks_real_random_dist8_sig15_15_noise_0.0025.mat",...
"0810_2\result_ADMM300iter_modes9_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.0025.mat",...
"0810_2\result_ADMM300iter_modes9_ort2_3ratio1stacks_real_random_dist8_sig15_15_noise_0.0025.mat"];
mydraw_compare_1(result_list,name_list)

% for what ? a suitable noise level - 0.675
mydraw_compare_1(result_list,name_list)

name_list = ["r=9,0","r=9,3"];
result_list = [
"0810_2\result_ADMM300iter_modes9_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat" ,...
"0810_2\result_ADMM300iter_modes9_ort2_3ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
];
mydraw_compare_1(result_list,name_list)

%% An important out come as 800 iter report
name_list = ["r=9","r=9,keep 6","r=6","r=6,without noise"];
result_list = [
"0811\result_ADMM800iter_modes9_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat" ,...
"0811\result_ADMM800iter_modes9_ort2_3ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
"0811\result_ADMM800iter_modes6_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
"0809_2\result_ADMM300iter_modes6_ort2_2ratio1stacks_real_random_dist8_sig15_15mat.mat"
];
mydraw_compare_1(result_list,name_list)

%0.25 no difference
%no notable difference
name_list = ["r=9,0","r=9,3","r=6,0","r=6,st"];
result_list = [
"0811noise\result_ADMM1000iter_modes9_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.03375.mat" ,...
"0811noise\result_ADMM1000iter_modes9_ort2_3ratio1stacks_real_random_dist8_sig15_15_noise_0.03375.mat",...
"0811noise\result_ADMM1000iter_modes6_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.03375.mat",...
"0809_2\result_ADMM300iter_modes6_ort2_2ratio1stacks_real_random_dist8_sig15_15mat.mat"
];
mydraw_compare_1(result_list,name_list)

%
name_list = ["1000","500","250","100","r=6,0","10"]; %problem
result_list = [
"0811noise\result_ADMM1000iter_modes6_ort2_2ratio1000beta0.05randdefaultstacks_real_random_dist8_sig15_15_noise_0.0675.mat" ,...
"0811noise\result_ADMM1000iter_modes6_ort2_2ratio500beta0.05randdefaultstacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
"0811noise\result_ADMM1000iter_modes6_ort2_2ratio250beta0.05randdefaultstacks_real_random_dist8_sig15_15_noise_0.0675.mat" ,...
"0811noise\result_ADMM1000iter_modes6_ort2_2ratio100beta0.05randdefaultstacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
"0811noise\result_ADMM800iter_modes6_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
"0811noise\result_ADMM1000iter_modes6_ort2_2ratio10beta0.05randdefaultstacks_real_random_dist8_sig15_15_noise_0.0675.mat"
];
mydraw_compare_1(result_list,name_list)

name_list = ["1000","500","r=9,0","r=9,3","r=6,0","r=6,st"]; %problem
result_list = [
"0811noise\result_ADMM1000iter_modes9_ort2_2ratio1000beta0.05randdefaultstacks_real_random_dist8_sig15_15_noise_0.0675.mat" ,...
"0811noise\result_ADMM1000iter_modes9_ort2_2ratio500beta0.05randdefaultstacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
"0811noise\result_ADMM800iter_modes9_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat" ,...
"0811noise\result_ADMM800iter_modes9_ort2_3ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
"0811noise\result_ADMM800iter_modes6_ort2_0ratio1stacks_real_random_dist8_sig15_15_noise_0.0675.mat",...
"0809_2\result_ADMM300iter_modes6_ort2_2ratio1stacks_real_random_dist8_sig15_15mat.mat"
];

%% compare data %report
data_list = [
     "stacks_real_random_dist8_sig15_15mat.mat",
     "stacks_random_dist8_blur.mat",
     "stacks_random_dist8_blur2.mat"
    ];

data_list = [
     "stacks_random_dist8_blur.mat",
     "stacks_real_random_dist8_blur_noise_0.000025.mat"
    ];


st_name = "stacks_real_random_st";
mydraw_compare_3(data_list,st_name);

%% new ADMM 
%r=1 r=3 r=6 
name_list = ["AP,r=12","ADMM,r=1","ADMM,r=3","ADMM r=9","ADMM,r=12"]
result_list = ["0806\APresult_ort0modes12stacks_random_dist8_blur.mat", ...,
    "0807\result_ADMM300iter_modes1_ort2_0stacks_random_dist8_blur.mat", ...,
    "0807\result_ADMM300iter_modes3_ort2_0stacks_random_dist8_blur.mat", ...,
    "0807\result_ADMM300iter_modes9_ort2_0stacks_random_dist8_blur.mat", ...,
    "0807\result_ADMM300iter_modes12_ort2_0stacks_random_dist8_blur.mat"]
mydraw_compare_1(result_list,name_list)

%% find optimal beta
% [5 1 5e-1 1e-1 5e-2 1e-2 5e-3 1e-3];
name_list = ["5", "1", "5e-1", "1e-1", "5e-2", "1e-2", "5e-3", "1e-3","old"];
result_list = ["0813\result_ADMM300iter_modes6_ort2_0ratio100beta5stacks_real_random_dist8_sig15_15mat.mat", ...,
    "0813\result_ADMM300iter_modes6_ort2_0ratio100beta1stacks_real_random_dist8_sig15_15mat.mat", ...,
    "0813\result_ADMM300iter_modes6_ort2_0ratio100beta0.5stacks_real_random_dist8_sig15_15mat.mat", ...,
    "0813\result_ADMM300iter_modes6_ort2_0ratio100beta0.1stacks_real_random_dist8_sig15_15mat.mat", ...,
    "0813\result_ADMM300iter_modes6_ort2_0ratio100beta0.05stacks_real_random_dist8_sig15_15mat.mat",...,
    "0813\result_ADMM300iter_modes6_ort2_0ratio100beta0.01stacks_real_random_dist8_sig15_15mat.mat",...,
    "0813\result_ADMM300iter_modes6_ort2_0ratio100beta0.005stacks_real_random_dist8_sig15_15mat.mat",...,
     "0813\result_ADMM300iter_modes6_ort2_0ratio100beta0.001stacks_real_random_dist8_sig15_15mat.mat",...,
     "0812\result_ADMM300iter_modes6_ort2_0ratio100stacks_real_random_dist8_sig15_15mat.mat"]
mydraw_compare_1(result_list,name_list)

name_list = [ "1e-2", "1e-3"];
result_list = [
    "0813\result_ADMM300iter_modes6_ort2_0ratio100beta0.01stacks_real_random_dist8_sig15_15mat.mat",...,
    "0813\result_ADMM300iter_modes6_ort2_0ratio100beta0.001stacks_real_random_dist8_sig15_15mat.mat"]
mydraw_compare_1(result_list,name_list)

%%
%find optimal beta2
name_list = ["1e-1", "5e-2", "1e-2", "5e-3", "1e-3"];
result_list = ["0908\result_ADMM500iter_modes12_ort2_0ratio100beta0.1stacks_regular_dist8_sig15_15.mat", ...,
    "0908\result_ADMM500iter_modes12_ort2_0ratio100beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0908\result_ADMM500iter_modes12_ort2_0ratio100beta0.01stacks_regular_dist8_sig15_15.mat", ...,
    "0908\result_ADMM500iter_modes12_ort2_0ratio100beta0.005stacks_regular_dist8_sig15_15.mat", ...,
    "0908\result_ADMM500iter_modes12_ort2_0ratio100beta0.001stacks_regular_dist8_sig15_15.mat"]
mydraw_compare_1(result_list,name_list)


name_list = ["100", "1000", "500", "250","10"];
result_list = ["0910\result_ADMM500iter_modes6_ort2_2ratio100beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio1000beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio500beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio250beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio10beta0.05stacks_regular_dist8_sig15_15.mat"]
mydraw_compare_1(result_list,name_list)

name_list = ["100", "1000", "500", "250","10"];
result_list = ["0910\result_ADMM500iter_modes3_ort2_2ratio100beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes3_ort2_2ratio1000beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes3_ort2_2ratio500beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes3_ort2_2ratio250beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes3_ort2_2ratio10beta0.05stacks_regular_dist8_sig15_15.mat"]
mydraw_compare_1(result_list,name_list)

name_list = ["100", "1000", "500", "250","10"];
result_list = ["0910\result_ADMM500iter_modes9_ort2_2ratio100beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes9_ort2_2ratio1000beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes9_ort2_2ratio500beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes9_ort2_2ratio250beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes9_ort2_2ratio10beta0.05stacks_regular_dist8_sig15_15.mat"]
mydraw_compare_1(result_list,name_list)

%% fig
name_list = ["1000", "500", "250","10","no ort","ort 1","ort 1 1"];
result_list = [
    "0910\result_ADMM500iter_modes6_ort2_2ratio1000beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio500beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio250beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio10beta0.05stacks_regular_dist8_sig15_15.mat",...,
    "0908modes\result_ADMM500iter_modes6_ort2_0ratio100beta0.05stacks_regular_dist8_sig15_15.mat",...,
    "0910\result_ADMM500iter_modes6_ort2_1ratio500beta0.05randdefaultstacks_regular_dist8_sig15_15.mat",...,
    "0928\result_ADMM300iter_modes6_ort2_1ratio1000beta0.05randdefaultstacks_regular_dist8_sig15_15.mat"
    ]
mydraw_compare_1(result_list,name_list)

name_list = ["1000", "500", "250","10","no ort"];
result_list = [
    "0910\result_ADMM501iter_modes6_ort2_2ratio1000beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio500beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio250beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0910\result_ADMM500iter_modes6_ort2_2ratio10beta0.05stacks_regular_dist8_sig15_15.mat",...,
    "0908modes\result_ADMM500iter_modes6_ort2_0ratio100beta0.05stacks_regular_dist8_sig15_15.mat"]
    
mydraw_compare_1(result_list,name_list,"model2_sig15_15.mat")


%% 0913 compare rand
%result_ADMM500iter_modes6_ort2_0ratio500beta0.05rand1stacks_regular_dist8_sig15_15.mat
name_list = ["6", "7", "8", "9","10"];
result_list = ["0914rand\result_ADMM500iter_modes6_ort2_0ratio500beta0.05rand1stacks_regular_dist8_sig15_15.mat", ...,
    "0914rand\result_ADMM500iter_modes6_ort2_0ratio500beta0.05rand2stacks_regular_dist8_sig15_15.mat", ...,
    "0914rand\result_ADMM500iter_modes6_ort2_0ratio500beta0.05rand3stacks_regular_dist8_sig15_15.mat", ...,
    "0914rand\result_ADMM500iter_modes6_ort2_0ratio500beta0.05rand4stacks_regular_dist8_sig15_15.mat", ...,
    "0914rand\result_ADMM500iter_modes6_ort2_0ratio500beta0.05rand5stacks_regular_dist8_sig15_15.mat"]
mydraw_compare_1(result_list,name_list)

name_list = ["6", "7", "8", "9","10"];
result_list = ["0914\result_ADMM500iter_modes6_ort2_2ratio500beta0.05rand6stacks_regular_dist8_sig15_15.mat", ...,
    "0914\result_ADMM500iter_modes6_ort2_2ratio500beta0.05rand7stacks_regular_dist8_sig15_15.mat", ...,
    "0914\result_ADMM500iter_modes6_ort2_2ratio500beta0.05rand8stacks_regular_dist8_sig15_15.mat", ...,
    "0914\result_ADMM500iter_modes6_ort2_2ratio500beta0.05rand9stacks_regular_dist8_sig15_15.mat", ...,
    "0914\result_ADMM500iter_modes6_ort2_2ratio500beta0.05rand10stacks_regular_dist8_sig15_15.mat"]
mydraw_compare_1(result_list,name_list)

%% compare modes _nes p1
c
%AP2000

%% compare approx p2
mydraw_compare_2("0908\result_ADMM500iter_modes3_ort2_0ratio100beta0.05stacks_regular_dist8_sig15_15.mat","model2_sig5_5.mat")
%% compare approx p3 4
result_list = ["0908modes\result_ADMM500iter_modes1_ort2_0ratio100beta0.05randdefaultstacks_regular_dist8_sig15_15.mat", ...,
    "0908modes\result_ADMM500iter_modes3_ort2_0ratio100beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0908modes\result_ADMM500iter_modes6_ort2_0ratio100beta0.05stacks_regular_dist8_sig15_15.mat", ...,
    "0908modes\result_ADMM500iter_modes9_ort2_0ratio100beta0.05stacks_regular_dist8_sig15_15.mat",...,
    "0908modes\result_ADMM500iter_modes12_ort2_0ratio100beta0.05stacks_regular_dist8_sig15_15.mat"]
mydraw_compare_4(result_list,"model2_sig15_15.mat")

%% solve problem
name_list = ["6", "7", "8", "9","10"];
result_list = ["0928\result_ADMM300iter_modes3_ort2_2ratio1000beta0.05rand0stacks_regular_dist8_sig15_15.mat", ...,
    "0928\result_ADMM300iter_modes3_ort2_2ratio1000beta0.05rand1stacks_regular_dist8_sig15_15.mat", ...,
    "0928\result_ADMM300iter_modes3_ort2_2ratio1000beta0.05rand2stacks_regular_dist8_sig15_15.mat", ...,
    "0928\result_ADMM300iter_modes3_ort2_2ratio1000beta0.05rand3stacks_regular_dist8_sig15_15.mat", ...,
    "0928\result_ADMM300iter_modes3_ort2_2ratio1000beta0.05rand4stacks_regular_dist8_sig15_15.mat"]
mydraw_compare_1(result_list,name_list)


%Interesting comparison
name_list = ["st","no st","real st"];
result_list = ["0928ini\result_ADMM300iter_modes3_ort2_0ratio1000beta0.05randdefaultstacks_regular_dist8_sig15_15.mat",...,
    "0928\result_ADMM300iter_modes3_ort2_1ratio1000beta0.05randdefaultstacks_regular_dist8_sig15_15.mat",...,
    "0928ini\result_ADMMnoblind300iter_modes3_ort2_0ratio1000beta0.05randdefaultstacks_regular_dist8_sig15_15.mat"]
mydraw_compare_1(result_list,name_list,"model2_sig15_15.mat")
mydraw_compare_4(result_list,"model2_sig15_15.mat")

%snrM and ort

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New mask
'data/stacks_regular_dist8_unregmask.mat'
load('data/stacks_regular_dist8_unregmask.mat')
ss = orthogonal(data.masks);
result.iter = 1;
result.masks = ss(:,:,1:20);
mydraw(result,1)

%%
load('data/stacks_regular_dist8_newmask.mat')
imshow(data.phobe,[])
result.masks = orthogonal(data.masks);
result.iter = 1;
result.masks = result.masks(:,:,1:4);
mydraw(result,1)    

phplot(result.masks(:,:,1))
axis equal
axis tight
title('\omega')

%%
[~, ~, rr]=make_grid(480,2);
r1 =5;
r2 = 19;
w = fftshift(rr>=r1 & rr <= r2);

aa = fftshift(ifft2(w));
aa2 = fftshift(ifft2(fftshift(w)));
imshowpair(aa,aa2,'montage');


