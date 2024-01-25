classdef Matrix
    properties
        dims
        dlist
        n
        r
        A
        B
        C
        U
        Up
    end
    
    methods
        function obj = Matrix(dlist, r)
            obj.dims = length(dlist);
            obj.dlist = dlist;
            obj.n = dlist(1);
            obj.r = r;
            obj.A = zeros(r, r);
            obj.B = zeros(r, r);
            obj.C = zeros(r, r);
            obj.U = cell(1, length(dlist));
            obj.Up = cell(1, length(dlist));
            for i = 1:length(dlist)
                obj.U{i} = zeros(dlist(i), r);
                obj.Up{i} = zeros(dlist(i), r);
            end
        end
        
        function Grad_F = get_proj_grad(obj, A_list, y)
            n = obj.dlist(1);
            sample_num = size(A_list, 2);
            y0 = zeros(sample_num, 1);
            for i = 1:sample_num
                a = reshape(A_list(:, i), n, 1);
                y0(i) = obj.A;
                for j = 1:obj.dims
                    y0(i) = y0(i) * (a' * obj.U{j});
                end
            end
            
            y0 = y0 - y;
            Grad_F = Matrix(obj.dlist, obj.r);
            for j = 1:obj.dims
                Grad_F.U{j} = reshape(obj.U{j}, obj.n, 1);
            end
            
            proj_list = zeros(obj.dims, sample_num);
            for i = 1:sample_num
                a = reshape(A_list(:, i), n, 1);
                for j = 1:obj.dims
                    proj_list(j, i) = a' * obj.U{j};
                end
            end
            
            Grad_F.A = prod(proj_list, 1) * y0;
            pos_list = true(1, obj.dims);
            for j = 1:obj.dims
                pos_list(j) = false;
                num = prod(proj_list(pos_list,:), 1); %careful! in Matlab a（[i]） will not get slices for first dim
                A_total = A_list * (y0 .* num');
                pos_list(j) = true;
                Grad_F.Up{j} = A_total - reshape(obj.U{j}, n, 1) * Grad_F.A;
            end
            
            return
        end

        function R_mat = retraction(obj, Grad_F, eta)
            R_mat = Matrix(obj.dlist, obj.r);
            Core = zeros(2 * ones(1, obj.dims));
            %index = ones(1, obj.dims); %careful when index Core([1,1]) is not a single number. should be Core(1,1)
            index = repmat({1}, 1, obj.dims);
            Core(index{:}) = obj.A - Grad_F.A * eta;
            for i = 1:obj.dims
                index{i} = 2;
                Core(index{:}) = -eta;
                index{i} = 1;
            end
        
            Qlist = cell(1, obj.dims);
            for j = 1:obj.dims
                UUp = [Grad_F.U{j}, Grad_F.Up{j}];
                [Qu, Ru] = qr(UUp, 0);
                Qlist{j} = Qu;
                Core = mul(Core, j, Ru);
            end
        
            B = Core;
            for j = 1:obj.dims
                mat = tensor_to_mat(Core, j);
                [u, s, vh] = svd(mat, 'econ');
                R_mat.U{j} = Qlist{j} * u(:, 1);
                B = mul(B, j, u(:, 1)');
            end
        
            R_mat.A = reshape(B, 1, 1);
            return
        end



    end
end