classdef Problem_CDP < Problem_PR % coded diffraction pattern problem
    properties
        % some data here
        L % number of diffraction patterns/ measurements
        N % number of pixels N by N figure
        D % list of random masks, size of L
        % p?
        m_id
        n_id
        % data about operator A: like maskD .
        % groundtruth x0
    end

    methods
        function obj = Problem_CDP(L,N,D,p)
            obj.L = L;
            obj.N = N;
            obj.D = D;
            obj.p = p;
            % constructor
            obj.m_id = [N,N,L];
            obj.n_id = [N,N]; 


        end

        function medium = aT(obj,R)  % aTR 
            % R single or multiphase(phase p) object. 
            % return: n by p matrix. rth row is a multiphase component
            % n is not necessary a scalar, may be N^2, or other shape
            % the same operation is taken on each phase
            % Replicate R to construct a obj.N x obj.N x obj.p x obj.L matrix
            medium = repmat(reshape(R, [obj.N, obj.N, obj.p, 1]), [1, 1, 1, obj.L]);
            medium = permute(medium,[1,2,4,3]);
            medium = myfft2(medium.*obj.D);
        end
        function result = a(obj,medium)  % sum a_r R^r 
            % medium single or multiphase(phase p) object.
            % medium [m_id, p] -> [n_id, p] 
            % return: n by p matrix. rth row is a multiphase component
            % n is not necessary a scalar, may be N^2, or other shape
            % the same operation is taken on each phase
            
            result = reshape((sum(myifft2(medium) .* conj(obj.D),3)),[n_id,p]); % sum over L mask, N by N by p

            %iFFT of lth mask, and sum over L mask?
        end

        % measurement, gradient, the same as Problem_PR
    end
end