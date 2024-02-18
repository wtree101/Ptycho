classdef Problem
   
    properties
        % some data here
        % data about operator A: like maskD .
        % groundtruth x0

    end

    methods
        function obj = Problem()
            % constructor
        end
        function measurement = A(obj,x) % forward operator
            % to be implemented
        end
        
        function gradient = Get_gradient(obj,x) %[A*(A(RR*)-b)]R
            % to be implemented
           
        end
    end
end

