n = 100;
r = 1;
dims = 3;
x0 = randn(n, r); % x0 is our target vector
x00 = x0/norm(x0);

m = 1000;
A_list = randn(n, m);
y = A(A_list, x0, dims);
X = Matrix([n,n,n], r);
n_iter = 1000;
eta = 1e-5;


err_list = zeros(1, n_iter);
for i = 1:n_iter
    grad_f = X.get_proj_grad(A_list, y);
    X = X.retraction(grad_f, eta);
    sign = x00' * X.U{1};
    err_list(i) = norm(x00 - sign * X.U{1});
end

plot(err_list);



function y0 = A(A_list, x0, dims)
    m = size(A_list, 2);
    y0 = zeros(m, 1);
    for i = 1:m
        a = A_list(:, i);
        y0(i) = (a.' * x0).^(dims);
    end
end

