modes= 15;
px=64; py=64;
load data2/stacks_real_regular_dist8_sig15_15.mat
amask = data.phobe; % normal amask
imshow(data.phobe,[])
mask_fre = fftshift(myfft2(ifftshift(data.phobe)));
%kernel_time = fspecial('average',5);
%kernel_time = fspecial('motion',20,45);
%kernel_time = padding(kernel_time,64);

kernel_fre = fftshift(myfft2( ifftshift(kernel_time)));

coherence = generate_toep(px,kernel_fre);
phobe = mask_fre(:)*mask_fre(:)';
density = coherence .* phobe;
[U,S,V] = svd(density);

q = U(:,1:modes)*S(1:modes,1:modes);
st_masks = reshape(q,[px py modes]);
st_masks = fftshift(myifft2(ifftshift(st_masks)));
result.masks = st_masks;
result.iter = 1;
mydraw(result,1);
