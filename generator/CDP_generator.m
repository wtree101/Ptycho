function CDP_generator(imageFlag,mask_num,problem)
    if ~exist('imageFlag','var')
        imageFlag=1;
    end
    scaleP=1;
    if imageFlag==1
        load GoldBalls2;a = A;clear A;
    end
    if imageFlag==2
        %a1=mat2gray(imread('cameraman.tif'));
        
        a0=mat2gray(imread('cameraman.tif'));
        a1=double(imresize(rgb2gray(imread('peppers.png')),[256,256]))/255;
        a1=a0.*exp(1j*2*pi*(a1)/2.*scaleP);
        
        a=min(1,abs(a1)).*sign(a1); %???camera man ??????????? ??????
        clear a0, a1;
    end
    
    if ~exist('sizeImg','var')
        sizeImg=size(a,1);
    end
    % sizeImg - to squarre
    a = imresize(real(a), [sizeImg sizeImg])+1j*imresize(imag(a), [sizeImg sizeImg]);
    R = a.reshape([size(a),1]); 
    Y = problem.A(R);

    data.sizeImg = sizeImg;
    data.image = a;
    data.nframes = mask_num;
    data.stacks = Y;
    %%
    % save
    save('data_CDP/CDP1.mat','data');
    

end
    %%