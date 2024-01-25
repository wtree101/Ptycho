classdef tensor_f
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
        function obj = tensor_f(dlist, r) %dlist = [m^2,M^2,M^2]
            obj.dims = length(dlist);
            obj.dlist = dlist;
            %obj.M = dlist(2);
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
        function obj_new = initialize(obj,ker,image)
            m = floor(sqrt(obj.dlist(1)));
            M = floor(sqrt(obj.dlist(2)));
            image = reshape(image,[M*M,1]);
            %ker = reshape(ker,[m*m,1]);
            obj_new = obj;
            obj_new.A = norm(image)^2 ;
            %obj_new.U{1} = ker/norm(ker);
            obj_new.U{1} = image/norm(image);
            obj_new.U{2} = conj(obj_new.U{1});

        end

        %function flj = forward_model(omega,position_j,)
        function phi = forward(obj,u,amask,cmapidx,is_conj)
            u_stack = u(cmapidx);
            
            if is_conj
                phi = conj(myfft2(conj(u_stack).*conj(amask)));
            else
                phi = myfft2(u_stack.* conj(amask));
            end
        end

        function fjq = myconv(obj,f,con_ker)
            newImage = flip(f,1);
            f_flip = flip(newImage,2); % f(-.)
            f_con = imfilter(f_flip,con_ker,'circular','conv'); %
            newImage = flip(f_con,1);
            fjq = flip(newImage,2); % f(-.)

        end
        
        function Grad_F = get_proj_grad(obj,cmapidx,y0,y,mask) % y:real_num y0:computed one 
            y0 = y0 - y;
            J = size(cmapidx,3);
            m = size(mask,1);
            M = floor(sqrt(obj.dlist(1)));
            Grad_F = tensor_f(obj.dlist, obj.r);
            for j = 1:obj.dims
                Grad_F.U{j} = obj.U{j};
            end

            y0_s = reshape(y0,[m,m,J]);
            %y0_shift_q = circshift(y0_s, [q(1), q(2)]);
%             newImage = flip(y0_s,1);
%             y0_s_flip = flip(newImage,2); % f(-.)

            %con_ker = reshape(Grad_F.U{1},[m,m]);
            V = reshape(Grad_F.U{1},[M,M]);
            W = reshape(Grad_F.U{2},[M,M]);
           
            FV = conj(obj.forward(V,mask,cmapidx,0));
            FW = conj(obj.forward(W,mask,cmapidx,1));
            
            

      
            product =  y0_s;
            Grad_F.A = sum(product(:));

           
            product = FW.*y0_s;
            ZV = FV;
            for j = 1:J
                ZV(:,:,j) = myifft2(product(:,:,j)) .* mask;
            end
            Qoverlap=@(frames) reshape(accumarray(cmapidx(:),frames(:),[M*M 1]),M,M);
            %when sum back, need to put patch back into bigger image:'
            Grad_F.Up{1} = reshape(Qoverlap(ZV),[M*M,1]) - obj.U{1} * Grad_F.A;
            
            %w
            product = FV.*y0_s;
            ZW = FW;
            for j = 1:J
                ZW(:,:,j) = conj(myifft2(conj(product(:,:,j)))) .* conj(mask);
            end
            Grad_F.Up{2} = reshape(Qoverlap(ZW),[M*M,1]) - obj.U{2} * Grad_F.A;
            return 

        end


        function R_mat = retraction(obj, Grad_F, eta)
            R_mat = tensor_f(obj.dlist, obj.r);
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
                %Qlist{j} = Qu(:,j);
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