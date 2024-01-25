function [density,masks,modes,ratio] = theoretical_model1(amask,porportion)
[px py] = size(amask);
kernel_fre = fspecial('gaussian',px,1);
%amask_time = myifft2(fftshift(myfft2(amask)));
%kernel_time = myifft2(kernel_fre);
%kernel_time = fftshift(kernel_time);
kernel_time = fftshift(myifft2(fftshift(kernel_fre)));
coherence = generate_toep(px,kernel_time);
phobe = amask(:)*amask(:)';
density = coherence .* phobe;
[U,S,V] = svd(density);
ratio =0; modes=0;
intensities = diag(S);
intensities  = intensities / sum(intensities(:));
while (ratio <= porportion)
    modes = modes + 1;
    ratio = ratio + intensities(modes);
end
q = U(:,1:modes)*S(1:modes,1:modes);
masks = reshape(q,[px py modes]);
end