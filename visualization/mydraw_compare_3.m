function mydraw_compare_3(data_list,st_name)
clf;
st = load(strcat("data\",st_name));
st = st.data; 
st = st.stacks;

num = size(data_list,1);
figure_ct = 1;

[px,py,nframes] = size(st);
medium = floor(nframes/3);

for i=1:num
    re{i} = load(strcat("data\",data_list(i)));
    re{i} = re{i}.data;
    re{i} = re{i}.stacks;
end

figure(figure_ct)
subplot(1,num+1,1);
imshow(st(:,:,medium),[]);
disp(['snr',num2str(snrCompt(st,st))]);
title("standard")
%figure_ct = figure_ct + 1;
%figure(figure_ct);
for i=1:num
    subplot(1,num+1,i+1);
    imshow(re{i}(:,:,medium),[])
    title(data_list(i),'Interpreter','none');
    %iter = re{i}.iter;
    disp([data_list(i),'snr',num2str(snrCompt((re{i}),(st)))]);
    %disp([name_list(i),'cor',num2str(coherence(re{i}.masks))]);
end

end