function mydraw_compare_1(result_list,name_list,standard)
%
%
% colormap(gray)


    
num = size(result_list,2);
%read data
if nargin==3
standard = load(strcat("standard\",standard));
standard = standard.standard;
st_masks = standard.st_masks;
end
for i=1:num
    re{i} = load(strcat("outcome\",result_list(i)));
    re{i} = re{i}.result;
    
end

% number

%compare figure
figure(1)
clf;
for i=1:num
    subplot(1,num,i);
    imshow(abs(re{i}.u),[])
    title(name_list(i));
    iter = re{i}.iter;
    disp([name_list(i),'R',num2str(re{i}.R(iter))]);
    disp([name_list(i),'cor',num2str(coherence(re{i}.masks))]);
    if nargin==3
    disp([name_list(i),'snrM', num2str(snrComptBlind(masks_to_rho(re{i}.masks),masks_to_rho(st_masks)))])
    end
    
end

figure(2)
clf;
hold on 
xlabel('Iteration number')
ylabel('R-factor')
%set(gca,'xscale','log')
%set(gca,'yscale','log')
%set(gca,'xscale','log')
% set(gca, 'XTick', [10 20 30 40 50 60 70 80])
% set(gca,'XTickLabel',{'1','10','25','50','100','150','200','1000'})

for i=1:num
    plot(re{i}.R,'DisplayName',name_list(i));
end
legend;
hold off

figure(3)
clf;
hold on 
xlabel('Iteration number')
ylabel('Snr')
for i=1:num
    plot(re{i}.snr,'DisplayName',name_list(i))
    
end
%set(gca,'xscale','log')
legend;
hold off

figure(4)
for i=1:num
    subplot(1,num,i);
    imshow(angle(re{i}.u),[])
    title(name_list(i));
end


figure(5)
clf;
hold on 
xlabel('Iteration number')
ylabel('Cor')
for i=1:num
    if isfield(re{i},'cor')
    plot(re{i}.cor,'DisplayName',name_list(i))
    end
end
%set(gca,'xscale','log')
legend;
hold off
% if isfield(result,'u')
%     u = result.u;
%     figure(1);
%     %subplot(2,modes,1)
%     imshow(abs(u));
%     title([num2str(result.iter),'_{th} relative error=',num2str(result.R(result.iter))]);
% end
% 
% if isfield(result,'masks')
%     masks = result.masks;
%     masks = orthogonal(masks);
%     [px,py,modes] = size(masks);
%     figure(2)
%     title('Time domain');
%     %colormap()
%     rows = floor((modes-1)/6) + 1;
%     %masks = fft2(masks);%fftshift(masks);
%     masks_fft = masks;
%     for k=1:modes
%         subplot(rows,6,k)
%         %imshow(abs(masks(:,:,k)),[]);
%         masks_fft(:,:,k) = myfft2(masks(:,:,k));
%         masks_fft(:,:,k) = fftshift(masks_fft(:,:,k));
%         % phplot(masks(:,:,k),1);
%         imshow(abs(masks(:,:,k)),[]);
%         title([num2str(k),'th mode'])
%         axis equal
%         axis tight
%     end
%     figure(3)
%     title('Fourier domain')
%     for k=1:modes
%         subplot(rows,6,k)
%         %imshow(abs(masks(:,:,k)),[]);
%         %     masks_fft(:,:,k) = myfft2(masks(:,:,k));
%         masks_fft(:,:,k) = fftshift(masks_fft(:,:,k));
%         % phplot(masks(:,:,k),1);
%         imshow(abs(masks_fft(:,:,k)),[]);
%         title([num2str(k),'th mode'])
%         axis equal
%         axis tight
%     end
%     
%     figure(4)
%     ss = reshape(masks,[px*py modes]);
%     q = diag(ss'*ss);
%     %q/sum(q)
%     s = sqrt(q)/sum(sqrt(q));
%     s = cumsum(s);
%     plot(s);
%     xlabel('Nth mode')
%     ylabel('Proportion(1~nth)')
%     title('Modes intensity')
% end
% 
% if isfield(result,'R')
%     figure(5)
%     plot(result.R);
% end
end