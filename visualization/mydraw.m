function mydraw(result,ort)
% colormap(gray)
close all
f_id = 1;
iter = result.iter;
   
    
if isfield(result,'u')
    u = result.u;
    figure(f_id);
    f_id = f_id + 1;
    %subplot(2,modes,1)
  %  imshow(abs(u),[]);
    phplot(u)
    title([num2str(result.iter),'_{th} relative error=',num2str(result.R(result.iter))]);
    figure(f_id);
    f_id = f_id + 1;
    %subplot(2,modes,1)
    imshow(angle(u),[]);
    title([num2str(result.iter),'_{th} relative error=',num2str(result.R(result.iter))]);
end

if isfield(result,'masks')
    masks = result.masks;
    if (ort==1)
    masks = orthogonal(masks);
    end
    [px,py,modes] = size(masks);
    figure(f_id);
    f_id = f_id + 1;
    %subplot(2,modes,1)
    
    title('Time domain');
    %colormap()
    rows = floor((modes-1)/12) + 1;
    %masks = fft2(masks);%fftshift(masks);
    masks_fft = masks;
    t = tiledlayout(rows,12);
    for k=1:modes
       % subplot(rows,12,k)
        %imshow(abs(masks(:,:,k)),[]);
        masks_fft(:,:,k) = myfft2(masks(:,:,k));
        masks_fft(:,:,k) = fftshift(masks_fft(:,:,k));
        % phplot(masks(:,:,k),1);
        nexttile
        imshow(abs(masks(:,:,k)),[],'border','tight','initialmagnification','fit');
        title([num2str(k),'th'])
        %axis equal
       % axis tight
        
    end
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    
    
    figure(f_id);
    f_id = f_id + 1;
    %subplot(2,modes,1)
    
    title('Time domain');
    %colormap()
    rows = floor((modes-1)/12) + 1;
    %masks = fft2(masks);%fftshift(masks);
    masks_fft = masks;
    for k=1:modes
        subplot(rows,12,k)
        %imshow(abs(masks(:,:,k)),[]);
        masks_fft(:,:,k) = myfft2(masks(:,:,k));
        masks_fft(:,:,k) = fftshift(masks_fft(:,:,k));
        % phplot(masks(:,:,k),1);
        %imshow(angle(masks(:,:,k)),[],'border','tight','initialmagnification','fit');
        phplot(masks(:,:,k))
        title([num2str(k),'th'])
        axis equal
        axis tight
    end
    
    
    figure(f_id);
    f_id = f_id + 1;
    %subplot(2,modes,1)
    
    title('Fourier domain')
    for k=1:modes
        subplot(rows,12,k)
        %imshow(abs(masks(:,:,k)),[]);
        %     masks_fft(:,:,k) = myfft2(masks(:,:,k));
        %masks_fft(:,:,k) = fftshift(masks_fft(:,:,k));
        % phplot(masks(:,:,k),1);
        imshow(abs(masks_fft(:,:,k)),[]);
        title([num2str(k),'th mode'])
        axis equal
        axis tight
    end
    
    figure(f_id);
    f_id = f_id + 1;
    ss = reshape(masks,[px*py, modes]);
    q = diag(ss'*ss);
    %q/sum(q)
    s = sqrt(q)/sum(sqrt(q));
   % s = cumsum(s);
    plot(s);
    xlabel('Nth mode')
    ylabel('Proportion(1~nth)')
    title('Modes intensity')
end

if isfield(result,'R')
    figure(f_id);
    f_id = f_id + 1;
    %subplot(2,modes,1)
    imshow(abs(u),[]);
    title([num2str(result.iter),'_{th} relative error=',num2str(result.R(result.iter))]);
    plot(result.R);
    disp(['R',num2str(result.R(iter))]);
end
if isfield(result,'snr')
   figure(f_id);
    f_id = f_id + 1;
    %subplot(2,modes,1)
    imshow(abs(u),[]);
    
    plot(result.snr);
    title([num2str(result.iter),'_{th} snr=',num2str(result.snr(result.iter))]);
end
if isfield(result,'cor')
    figure(f_id);
    f_id = f_id + 1;
    %subplot(2,modes,1)
    imshow(abs(u),[]);
    
    plot(result.cor);
    title([num2str(result.iter),'_{th} cor=',num2str(result.cor(result.iter))]);
    disp(['cor',num2str(coherence(result.masks))]);
end

end