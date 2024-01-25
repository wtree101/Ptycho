function [density,st_masks,modes,ratio,intensities] = theoretical_model1_revised(amask,porportion)
[px py] = size(amask);
name = ['standard\model1_disk5_pro',num2str(porportion),'.mat'];
kernel_fre = fspecial('guassian',64,1);
%kernel_fre = padding(fspecial('disk',5),px);
%amask_time = myifft2(fftshift(myfft2(amask)));
kernel_time = myifft2(kernel_fre);
kernel_time = fftshift(kernel_time);
coherence = generate_toep(px,kernel_time);
phobe = amask(:)*amask(:)';
density = coherence .* phobe;
[U,S,V] = svd(density);
intensities = diag(S);
intensities  = intensities / sum(intensities(:));
if (porportion<1)
    ratio =0; modes=0;
    
    while (ratio <= porportion && modes <= px*px)
        modes = modes + 1;
        ratio = ratio + intensities(modes);
    end
else
    modes=px*px;
    ratio = 1;
end
q = U(:,1:modes)*S(1:modes,1:modes);
st_masks = reshape(q,[px py modes]);
%st_density = density;
%save([name,'density'],'st_density');
%save(name,'st_masks');
end