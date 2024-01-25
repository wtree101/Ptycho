function mydraw_compare_4(result_list,standard)
%
%
% colormap(gray)
%num = size(result_list,2);
%read data
close all

num = size(result_list,2);
%read data

for i=1:num
    re{i} = load(strcat("outcome\",result_list(i)));
    re{i} = re{i}.result;
    
end

standard = load(strcat("standard\",standard));
standard = standard.standard;
st_masks = standard.st_masks;
[px,py,st_modes] = size(st_masks);

standard_density = standard.st_q;
%[u_st,s_st,v_st] = svd(standard_density,'econ');
st_s = standard.st_s;
st_u = standard.st_u;
%compare figure



f_id = 0;


%%

f_id = f_id + 1;
figure(f_id)
%ss_st = reshape(standard,[px*py modes]);
%q = diag(ss_st'*ss_st);
%q/sum(q)

s = (st_s.^2)/sum(st_s.^2);
s = cumsum(s);
plot(s(1:20));
xlabel('Nth mode')
ylabel('Proportion(1~nth)')
title('Standard modes intensity(accumulative)')
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

s_2 = st_s/st_s(1);
plot(s_2(1:20));
xlabel('Nth mode')
ylabel('s_i / s_1')
title('Standard mode intensity')

%%
f_id = f_id + 1;
figure(f_id)
%ss_st = reshape(standard,[px*py modes]);
%q = diag(ss_st'*ss_st);
%q/sum(q)

hold on
errM_opt = sqrt(1 - s(1:12));
errM = zeros(1,num);
modes_list = zeros(1,num);
plot(errM_opt,'DisplayName','err_M^*');

for i=1:num
    result = re{i};
    masks = result.masks;
%     masks = orthogonal(masks,1);
    [px,py,modes] = size(masks);
    ss = reshape(masks,[px*py modes]);
    density = ss*ss';
    scale=sum(standard_density(:).*conj(density(:)))/norm(density,'fro')^2;
    density = density * scale;
    density_err = norm(density - standard_density,'fro')/norm(standard_density,'fro');
    errM(i) = density_err;
    modes_list(i) = modes;
end

plot(modes_list,errM,'DisplayName','err_M');

legend;

xlabel('N modes')
ylabel('Approximation error err_M')
title('Approximate the standard density matrix')
%
% if isfield(result,'R')
%     figure(5)
%     plot(result.R);
% end

%%



end

