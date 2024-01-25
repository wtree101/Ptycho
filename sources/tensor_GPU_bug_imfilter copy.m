classdef tensor_GPU
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
        type
    end
    
    methods
        function obj = tensor_GPU(dlist, r, type) %dlist = [m^2,M^2,M^2]
            if nargin < 3
                type = 'CPU';
            end
            obj.dims = length(dlist);
            obj.dlist = dlist;
            obj.type = type;
            if type=="CPU"
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
            else
                %GPU version
                obj.r = gpuArray(r);
                obj.A = gpuArray.zeros(r, r);
                obj.B = gpuArray.zeros(r, r);
                obj.C = gpuArray.zeros(r, r);
                obj.U = cell(1, length(dlist));
                obj.Up = cell(1, length(dlist));
                for i = 1:length(dlist)
                    obj.U{i} = gpuArray.zeros(dlist(i), r);
                    obj.Up{i} = gpuArray.zeros(dlist(i), r);
                end
            end

            
        end
        function obj_new = initialize(obj,ker,image)
            m = floor(sqrt(obj.dlist(1)));
            M = floor(sqrt(obj.dlist(2)));
            image = reshape(image,[M*M,1]);
            ker = reshape(ker,[m*m,1]);
            obj_new = obj;
            obj_new.A = norm(image)^2 * norm(ker);
            obj_new.U{1} = ker/norm(ker);
            obj_new.U{2} = image/norm(image);
            obj_new.U{3} = conj(obj_new.U{2});

            if obj.type=="GPU"
                obj_new.A = gpuArray(norm(image)^2 * norm(ker));
                obj_new.U{1} = gpuArray(ker/norm(ker));
                obj_new.U{2} = gpuArray(image/norm(image));
                obj_new.U{3} = gpuArray(conj(obj_new.U{2}));
            end

        end

        %function flj = forward_model(omega,position_j,)
        function phi = forward(obj,u,amask,cmapidx,is_conj)
            % assume all var are on GPU or CPU
            u_stack = u(cmapidx);
            
            if is_conj
                phi = conj(myfft2(conj(u_stack).* amask)); %QU ?
            else
                phi = myfft2(u_stack.* amask); %QU ?
            end
        end

        function fjq = myconv(obj,f,con_ker)
            newImage = flip(f,1);
            f_flip = flip(newImage,2); % f(-.)
            f_con = imfilter(f_flip,con_ker,'circular','conv'); %
            newImage = flip(f_con,1);
            fjq = flip(newImage,2); % f(-.)

        end

        function Grad_F = get_proj_grad(obj,cmapidx,FV,FW,FVFW,y_now,y,mask,Qoverlap) % y:real_num y0:computed one 
            %assert all inputs are on CPU or GPU
            y0 = y_now - y;
            J = size(cmapidx,3);
            m = floor(sqrt(obj.dlist(1)));
            M = floor(sqrt(obj.dlist(2)));
            if obj.type=="CPU"
                Grad_F = tensor_GPU(obj.dlist, obj.r);
            else
                Grad_F = tensor_GPU(obj.dlist, obj.r,"GPU");
            end

            for j = 1:obj.dims
                Grad_F.U{j} = obj.U{j};
            end

            y0_s = reshape(y0,[m,m,J]);
            %y0_shift_q = circshift(y0_s, [q(1), q(2)]);
%             newImage = flip(y0_s,1);
%             y0_s_flip = flip(newImage,2); % f(-.)

            con_ker = reshape(Grad_F.U{1},[m,m]);
            V = reshape(Grad_F.U{2},[M,M]);
            W = reshape(Grad_F.U{3},[M,M]);
           
            %FV = conj(obj.forward(V,mask,cmapidx,0));
            %FW = conj(obj.forward(W,mask,cmapidx,1));
            
            
%             con_f_pre = imfilter(y0_s_flip,con_ker,'circular','conv'); %
%             newImage = flip(con_f_pre,1);

            %con_f = flip(newImage,2); % .. (-q)
%             con_f = y0_s;
%             parfor j = 1:J
%                 %y0_I(:,:,j) = imfilter(y0_I(:,:,j),product(:,:,j),'circular','conv' );
%                 con_f(:,:,j) = obj.myconv(conj(y0_s(:,:,j)),con_ker);
%             end
%             con_f = conj(con_f);

            %con_f = conj(obj.myconv(conj(y0_s),con_ker));   %不知道conv这些对前两个维度操作会怎么样。保险还是分开写
            con_f = conj(imfilter(conj(y0_s),con_ker,'circular' ));
            %con_f2 = obj.myconv_vectorized(conj(y0_s), con_ker);
      
            product = FVFW.*con_f;
            Grad_F.A = sum(product(:));


            % con_ker p
            product = FVFW;
            y0_I = y0_s;
            
            parfor j = 1:J
                %tic
                y0_I(:,:,j) = imfilter(y0_s(:,:,j),product(:,:,j),'circular');
                %toc
                %y0_I(:,:,j) = obj.myconv(y0_s(:,:,j),product(:,:,j));
            end
            
          %  y0_I_2 = obj.myconv(y0_s,product);


            

            Grad_F.Up{1} = reshape(sum(y0_I,3),[m*m,1]) - obj.U{1} * Grad_F.A;

            % v
            product = con_f.*FW;
%             ZV = FV;
%             parfor j = 1:J
%                 ZV(:,:,j) = myifft2(product(:,:,j)) .* conj(mask); %QU ? corespond to ADMM
%             end

            ZV = myifft2(product) .* conj(mask);
            %Qoverlap=@(frames) reshape(accumarray(cmapidx(:),frames(:),[M*M 1]),M,M);
            %when sum back, need to put patch back into bigger image:'
            Grad_F.Up{2} = reshape(Qoverlap(ZV),[M*M,1]) - obj.U{2} * Grad_F.A;
            
            %w
            product = FV.* con_f;
%             ZW = FW;
%             parfor j = 1:J
%                 ZW(:,:,j) = conj(myifft2(conj(product(:,:,j)))) .* mask ;
%             end

            ZW = conj(myifft2(conj(product))) .* mask;
            Grad_F.Up{3} = reshape(Qoverlap(ZW),[M*M,1]) - obj.U{3} * Grad_F.A;
            

            return 

        end


        function R_mat = retraction(obj, Grad_F, eta)
            % assert Grad_F and obj are on the same place
            if obj.type=="CPU"
                R_mat = tensor_GPU(obj.dlist, obj.r);
                Core = zeros(2 * ones(1, obj.dims));
            else
                R_mat = tensor_GPU(obj.dlist, obj.r,"GPU");
                Core = gpuArray.zeros(2 * ones(1, obj.dims), 'single');
            end
            
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
