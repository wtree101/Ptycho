function [masks_new] = orthogonal_reduce(masks,mode_keep)
[px,py,modes] = size(masks);

ss = reshape(masks,[px*py modes]);
%azzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz [u,s,~] = svd(ss*ss');
% q = u(:,1:modes)*s(1:modes,1:modes);
% [V,D] = eig(ss*ss');
% [d,ind] = sort(diag(D),'descend');
% ind = ind(1:modes);
% Ds = D(ind,ind);
% Vs = V(:,ind);
% q = Vs*sqrt(Ds);
[U,S,V] = svd(ss,'econ');
S = diag(S);
S(mode_keep+1:end)=0;
S = diag(S);
q = U*S*V';
masks_new = reshape(q,[px py modes]);


end