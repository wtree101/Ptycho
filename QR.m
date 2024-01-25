%using PyPlot
A = rand(10,10); 
[q,r] = qr(A);
A = q' * diag([512,256,128,64,32,16,8,4,2,1]) * q;
for k = 1:100000
    mu = A(mod(k,10)+1,mod(k,10)+1);
    [Q,R] = qr(A-mu*eye(10));
    A = R*Q + mu*eye(10);
    A
end
spy(abs(A)>1e-4)