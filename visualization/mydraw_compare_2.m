function mydraw_compare_2(name,standard)
%
%
% colormap(gray)
%num = size(result_list,2);
%read data
close all

result = load(strcat("outcome\",name));
result = result.result;
standard = load(strcat("standard\",standard));
standard = standard.standard;
st_masks = standard.st_masks;
[px,py,st_modes] = size(st_masks);

standard_density = standard.st_q;
%[u_st,s_st,v_st] = svd(standard_density,'econ');
st_s = standard.st_s;
st_u = standard.st_u;
%compare figure

masks = result.masks;
masks = orthogonal(masks);
[px,py,modes] = size(masks);

f_id = 1;

figure(f_id)
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
    imshow(abs(masks(:,:,k)),[]);
    title([num2str(k),'th mode'])
    axis equal
    axis tight
end

f_id = f_id + 1;
figure(f_id)
%standard = orthogonal(standard);
for k=1:modes
    subplot(rows,12,k)
    %imshow(abs(masks(:,:,k)),[]);
    %     masks_fft(:,:,k) = myfft2(masks(:,:,k))
    % phplot(masks(:,:,k),1);
    imshow(abs(st_masks(:,:,k)),[]);
    title([num2str(k),'th mode'])
    axis equal
    axis tight
end

% for k=1:modes
%     subplot(rows,6,k)
%     %imshow(abs(masks(:,:,k)),[]);
%     %     masks_fft(:,:,k) = myfft2(masks(:,:,k));
%     masks_fft(:,:,k) = fftshift(masks_fft(:,:,k));
%     % phplot(masks(:,:,k),1);
%     imshow(abs(masks_fft(:,:,k)),[]);
%     title([num2str(k),'th mode'])
%     axis equal
%     axis tight
% end


f_id = f_id + 1;
figure(f_id)
ss = reshape(masks,[px*py modes]);
q = diag(ss'*ss);
%q/sum(q)
s = sqrt(q)/sum(sqrt(q));
s = cumsum(s);
plot(s);
xlabel('Nth mode')
ylabel('Proportion(1~nth)')
title('Modes intensity')



f_id = f_id + 1;
figure(f_id)
%ss_st = reshape(standard,[px*py modes]);
%q = diag(ss_st'*ss_st);
%q/sum(q)

s = (st_s.^2)/sum(st_s.^2);
s = cumsum(s);
plot(s);
xlabel('Nth mode')
ylabel('Proportion(1~nth)')
title('Modes intensity st')
%
% if isfield(result,'R')
%     figure(5)
%     plot(result.R);
% end
%%
f_id = f_id + 1;
figure(f_id)
%ss_st = reshape(standard,[px*py modes]);
%q = diag(ss_st'*ss_st);
%q/sum(q)

plot(st_s);
xlabel('Nth mode')
ylabel('Proportion(1~nth)')
title('Modes intensity st')

%% density comparison
f_id = f_id + 1;
figure(f_id)

subplot(1,3,1)
imshow(abs(standard_density),[]);
%title(['standard-approx vs density: ',num2str(density_err),'  singular value',num2str(sqrt(1-s(modes)))])
title('standard density')

subplot(1,3,2)
st_S = diag(st_s);
%q = st_u(:,1:modes)*sqrt(st_S(1:modes,1:modes);
standard_density_approx = st_u(:,1:modes)*st_S(1:modes,1:modes)*st_u(:,1:modes)';
density_err = norm(standard_density - standard_density_approx,'fro')/norm(standard_density,'fro');
imshow(abs(standard_density_approx),[]);
title(['standard vs standard-approx: ',num2str(density_err),'  singular value',num2str(sqrt(1-s(modes)))])

subplot(1,3,3)
density = ss * ss';
scale=sum(standard_density(:).*conj(density(:)))/norm(density,'fro')^2;
density = density * scale;
density_err = norm(density - standard_density_approx,'fro')/norm(standard_density,'fro');
imshow(abs(density),[]);
title(['standard-approx vs result: ',num2str(density_err)])

density_err = norm(standard_density - density,'fro')/norm(standard_density,'fro');
sgtitle(['standard vs result: ',num2str(density_err)])

%%



end

