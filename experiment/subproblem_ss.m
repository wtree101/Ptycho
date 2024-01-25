function [rho,dist] = find_min_rho(B,iter)
[n,r] = size(B);
alpha = ones(1,r);

for i = 1:iter
   [U S V] = svd( B * diag(alpha)' );
   % omega update
   omega = U * V';
   % alpha update
   alpha = sum(B,1) ./ sum(omega,1);
   
end
rho = omega * diag(alpha);
dist = form(rho - B,'fro');

end