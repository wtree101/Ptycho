modes= 15;
px=64; py=64;
load 'data/stacks_regular_dist8_blur1_new.mat'
amask = data.phobe; % normal amask
imshow(data.phobe,[])
%mask_fre = fftshift(myfft2(ifftshift(data.phobe)));
kernel_fre = fspecial('gaussian',64,1);

kernel_time = fftshift(myifft2( ifftshift(kernel_fre)));

coherence = generate_toep(px,kernel_time);
phobe = amask(:)*amask(:)';
density = coherence .* phobe;
[U,S,V] = svd(density);

q = U(:,1:modes)*S(1:modes,1:modes);
st_masks = reshape(q,[px py modes]);

result.masks = st_masks;
result.iter = 0;
mydraw(result,1);
