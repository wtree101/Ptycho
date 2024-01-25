function [rho,dist] = find_min_rho(B,iter)
[n,r] = size(B);
alpha = ones(1,r);

for i = 1:iter
   [U S V] = svd( B * diag(alpha)' ,'econ');
   % omega update
   omega = U * V';
   % alpha update
   %alpha = sum(B,1) ./ sum(omega,1);
   alpha = sum( real(conj(omega).* B),1);  %./ sum( abs(omega).^2,1); = 1 
   
end
rho = omega * diag(alpha);
dist = norm(rho - B,'fro');

end