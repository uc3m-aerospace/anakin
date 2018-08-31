%{
inertia: class to define inertia tensors. Inherits from tensor

I0 = anakin.inertia();  % return default object
I  = anakin.inertia(I|T|m,<a|c>,<S1>);

where: 
- <> denotes optional arguments
- | denotes alternative arguments
- () groups argument options
- I0 is the default inertia object(null tensor)
- I is a inertia object 
- T is a tensor   
- m is a matrix
- a is a vector
- c is an array with the three vector components
- S1 is a frame. If given, all previous input as relative to that frame
 
METHODS:
* matrix: returns the matrix of the tensor in a chosen frame
* eigs: returns eigenvectors and eigenvalues 
* steiner: displaces the point about which the tensor....????
* dt: returns the time derivative of a tensor wrt a chosen frame
  (symbolic variables must be used)  

AUTHOR: Mario Merino <mario.merino@uc3m.es>
%}
classdef inertia < anakin.tensor  
    methods % creation
        function I = inertia(varargin) % constructor
            for i = 1:length(varargin)
               if isa(varargin{i},'sym')
                   varargin{i} = formula(varargin{i}); % enforce formula to allow indexing
               end
            end            
            switch nargin
                case 0 % no arguments
                    return;
                case 1  
                    I.m = anakin.inertia(varargin{1},anakin.frame).m;
                case 2 
                    if isa(varargin{end},'anakin.frame') % last is frame
                        if isa(varargin{1},'anakin.tensor') % tensor, frame
                            I.m = varargin{2}.matrix * varargin{1}.m * varargin{2}.matrix';
                        else % matrix, basis
                            I.m = varargin{2}.matrix * varargin{1} * varargin{2}.matrix';
                        end
                    else
                        I.m = anakin.inertia(varargin{1},varargin{2},anakin.frame).m;
                    end    
                case 3
                    I.m = varargin{3}.matrix * varargin{1}.m * varargin{3}.matrix';
                    I = I.steiner(anakin.vector(-anakin.vector(varargin{2}).c-varargin{3}.c));
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in tensor');
            end       
        end
    end
    methods % overloads
        function value = eq(I1,I2) % overload ==
            if isa(I1.m,'sym') || isa(I1.m,'sym') % symbolic inputs
                value = isAlways(I1.m==I2.m,'Unknown','false'); % In case of doubt, false
            else % numeric input
                value = (abs(I1.m - I2.m) < 10*eps(I1.m)+10*eps(I2.m)); 
            end
            value = all(value(:));
        end
        function value = ne(I1,I2) % overload ~=
            value = ~eq(I1,I2);
        end
        function I1 = plus(I1,I2) % overloaded + operator
            I1.m = I1.m + I2.m; 
        end
        function I1 = minus(I1,I2) % overloaded - operator
            I1.m = I1.m - I2.m; 
        end 
        function I1 = uplus(I1) % overloaded + operator (unitary)
            % pass
        end
        function I1 = uminus(I1) % overloaded - operator (unitary)
            I1.m = -I1.m; 
        end 
        function I1 = times(I1,I2) % overloaded .* (multiplication by scalar)
            if isa(I1,'anakin.tensor') % then b is not ten
                I1.m = I1.m .* I2; 
            else % a is not vector
                I2.m = I1 .* I2.m; 
                I1 = I2;
            end                
        end
        function I1 = mtimes(I1,I2) % overloaded * (multiplication by scalar or matrix)
            if isa(I1,'anakin.tensor') % then b is not vector
                I1.m = I1.m * I2; 
            else % a is not vector
                I2.m = I1 * I2.m; 
                I1 = I2;
            end
        end 
        function I1 = rdivide(I1,x) % overloaded ./ (division by scalar)
            I1.m = I1.m ./ x; 
        end
        function I1 = mrdivide(I1,x) % overloaded / (division by scalar or matrix)
            I1.m = I1.m / x; 
        end
        function I1 = ldivide(x,I1) % overloaded .\ (division by scalar)
            I1.m = x .\ I1.m; 
        end
        function I1 = mldivide(x,I1) % overloaded \ (division by scalar or matrix)
            I1.m = x \ I1.m; 
        end         
        function disp(I) % display
            disp('Inertia tensor with canonical components')
            disp('with respect to the canonical origin:')
            disp(I.m)
        end
    end 
    methods % functionality
        function matrix = matrix(I,B1) % matrix of the tensor in B1
            if ~exist('B1','var')
                B1 = anakin.basis; % default frame
            end 
            matrix = B1.m' * I.m * B1.m;
            if isa(matrix,'sym')
                matrix = formula(simplify(matrix));
            end
        end
        function [l1,l2,l3,v1,v2,v3] = eigs(I) % like Matlab eigs function
           [V,D] = eigs(I.m);
           l1 = D(1,1);
           l2 = D(2,2);
           l3 = D(3,3);           
           v1 = anakin.vector(V(:,1));
           v2 = anakin.vector(V(:,2));
           v3 = anakin.vector(V(:,3));
        end
    end
    methods % symbolic
        function dI = dt(I,B) % time derivative with respect to basis B. Requires sym vector that utlimately depends on a single variable t
            if ~exist('B','var')
                B = anakin.basis; % canonical vector basis
            end
            dI = anakin.tensor(diff(sym(I.matrix(B)),1),B);
        end
        function I_ = subs(I,variables,values) % particularize symbolic vector
            I_ = I;
            I_.c = double(subs(I.m,variables,values));
        end
    end
    methods % logical tests
        function ishermitian = ishermitian(I) % tensor is hermitian (A' * A = eye(3))
            if isa(a.c,'sym') % symbolic inputs
                ishermitian = isAlways(I.m' * I.m == eye(3),'Unknown','false'); % In case of doubt, false
            else % numeric input            
                ishermitian = (abs(I.m' * I.m - eye(3))<eps(max(abs(I.m)))); 
            end 
            ishermitian = all(ishermitian);
        end 
    end
end








