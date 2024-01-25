load("outcome\result_ADMM_200iter_modes15_ort0stacks_random_dist8_blur.mat")
masks = result.masks;
main = masks(:,:,1);
mydraw(result)
imshow(abs(diff(main,1,2)),[])
%imshow(abs(diff(main,1,1)),[])
%mydraw(result)