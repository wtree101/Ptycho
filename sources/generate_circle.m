function fmask = generate_circle(px,r)
sizeProb = px;
scale_r=floor(sqrt(35/3*5/3*sizeProb/64.)* r * 2);
ker0 = fspecial('disk',scale_r);
fmask = padding(ker0,px);
end
