load 'data/stacks_regular_dist8_blur1_new.mat'
amask = data.phobe;
px=64;
py=64;
modes=2;
masks = amask(:,:,ones(2,1));
masks(:,:,1) = circshift(amask,[4,4]);
masks(:,:,2) = circshift(amask,[-4,-4]);
% masks(:,:,4) = circshift(amask,[4,0]);
% 
% masks(:,:,5) = circshift(amask,[-4,0]);
% masks(:,:,6) = circshift(amask,[4,4]);
% masks(:,:,7) = circshift(amask,[4,-4]);
% masks(:,:,8) = circshift(amask,[-4,4]);
% masks(:,:,9) = circshift(amask,[-4,-4]);

dd = reshape(masks,[px*py modes]);

result.iter=1;
result.masks=masks;
mydraw(result,1);