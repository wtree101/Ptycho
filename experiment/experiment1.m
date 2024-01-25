gu = fspecial('gaussian',64,1);
%amask_time = myifft2(fftshift(myfft2(amask)));
ss = myifft2(gu);
ss = fftshift(ss);
coherence = generate_toep(64,ss);
phobe = 
q = U(:,1:12)*S(1:12,1:12);