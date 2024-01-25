function rho = masks_to_rho(masks)
[px,py,modes] = size(masks);
ss = reshape(masks,[px*py modes]);
rho = ss * ss';
end