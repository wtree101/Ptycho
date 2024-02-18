classdef Problem_PR < Problem
    % class for the problem of phase retrieval
   
    properties
        % some data here
        L % number of diffraction patterns/ measurements
        N % number of pixels N by N figure
        D % list of random masks, size of L
        % p?
        m_id
    

    end
    % data about operator A: like maskD .
    % groundtruth x0

    methods
        function obj = Problem_PR()
            % constructor
            obj.p = 1;
        end
        function medium = aT(obj,R)  % aTR 
            % R single or multiphase(phase p) object. 
            % return: n by p matrix. rth row is a multiphase component
            % n is not necessary a scalar, may be N^2, or other shape
            % the same operation is taken on each phase
        end
        function result = a(obj,medium)  % sum a_r R^r 
            % R single or multiphase(phase p) object. 
            % return: n by p matrix. rth row is a multiphase component
            % n is not necessary a scalar, may be N^2, or other shape
            % the same operation is taken on each phase
        end
        function measurement = A(obj,R)  %A(RRT)
            % forward operator
            % for PR problem, measurement comes from square of abs of medium
            % for multi-phase, measurement is the sum of square of abs of each phase
            %idx = repmat({':'}, 1, ndims(obj.m_id));
            medium = obj.aT(R);
            % Combine the measurement of each phase using matrix operations
            % for single phase, make sure last dim is 1(should keep it in the medium)
            measurement = sum(abs(medium).^2, ndims(medium)); %last dim, squeeze is not needed
            % to be checked: ?
            % usage of idx, ndims
        end
        
        function gradient = Get_gradient(obj,R) %[A*(A(RR*)-b)]R
            % sum (y_r a^r) (medium^r)  
            % for multi-phase, each phase takes the same operation.
            medium = obj.aT(R);
            y0 = obj.A(R); 
            % AT(y0) R 
            medium = medium .* y0; % act on each phase ?
            gradient = a(medium);
        end
    end
end