function pre_standard(input,type,modes_keep)
%
%
% colormap(gray)
%num = size(result_list,2);
%read data

if (type=="masks")
    [px,py,st_modes] = size(input);
    q_st = reshape(input, [px*py st_modes]);
    standard_density = q_st * q_st';
else
    standard_density = input;
    st_modes = size(input,1);
    px = floor(sqrt(st_modes));
    
end

if nargin==2
    modes_keep = st_modes;
end

[u_st,s_st,v_st] = svd(standard_density,'econ');
%s_st = sqrt(s_st);
standard.st_q = standard_density;
standard.st_s = diag(s_st(1:modes_keep,1:modes_keep));
standard.st_u = u_st(:,1:modes_keep);
%q = u_st(:,1:modes_keep)*s_st(1:modes_keep,1:modes_keep);
q =  u_st(:,1:modes_keep)*sqrt(s_st(1:modes_keep,1:modes_keep));
st_masks = reshape(q,[px py modes_keep]);
standard.st_masks = st_masks;

%compare figure
save('model\standard\model1_gu1','standard');

end