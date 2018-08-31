%{
tensor: class to define 3x3-tensors.

T0 = anakin.tensor();  % returns default object
T  = anakin.tensor(T|m|((a|c),(a|c)),<B1>);

where: 
- <> denotes optional arguments
- | denotes alternative arguments
- () groups argument options
- T0 is the default tensor (null tensor)
- T is a tensor  
- m is a matrix
- a is a vector
- c is an array with the three vector components
- B1 is a basis. If given, all previous input as relative to that basis
 
METHODS:
* matrix: returns the matrix of the tensor in a chosen basis
* eigs: returns eigenvectors and eigenvalues
* isunitary: checks if tensor is unitary (Hermitian)
* dt: returns the time derivative of a tensor wrt a chosen basis
  (symbolic variables must be used)  
* subs: takes values of the symbolic unknowns and returns a vector with
  purely numeric coordinates (symbolic variables must be used)   

AUTHOR: Mario Merino <mario.merino@uc3m.es>
%}
classdef tensor
    properties (Hidden = true, Access = protected) 
        m = [0,0,0;0,0,0;0,0,0]; % tensor components in the canonical vector basis B0
    end 
    methods % creation
        function T = tensor(varargin) % constructor
            for i = 1:length(varargin)
               if isa(varargin{i},'sym')
                   varargin{i} = formula(varargin{i}); % enforce formula to allow indexing
               end
            end            
            switch nargin
                case 0 % no arguments
                    return;
                case 1  
                    T.m = anakin.tensor(varargin{1},anakin.basis).m;
                case 2 
                    if isa(varargin{end},'anakin.basis') % last is basis
                        if isa(varargin{1},'anakin.tensor') % tensor, basis
                            T.m = varargin{2}.matrix * varargin{1}.m * varargin{2}.matrix';
                        else % matrix, basis
                            T.m = varargin{2}.matrix * varargin{1} * varargin{2}.matrix';
                        end
                    else
                        T.m = anakin.tensor(varargin{1},varargin{2},anakin.basis).m;
                    end   
                case 3 % vector, vector, basis
                    T.m = tensorproduct(anakin.vector(varargin{1},varargin{3}),anakin.vector(varargin{2},varargin{3}));
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in tensor');
            end       
        end
        function T = set.m(T,value) % on setting c
            T.m = reshape(value,3,3); % Force column 
            if isa(T.m,'sym') % symbolic input
                T.m = formula(simplify(T.m)); % simplify and force sym rather than symfun to allow indexing into c
            end
        end
    end
    methods % overloads
        function value = eq(T1,T2) % overload ==
            if isa(T1.m,'sym') || isa(T1.m,'sym') % symbolic inputs
                value = isAlways(T1.m==T2.m,'Unknown','false'); % In case of doubt, false
            else % numeric input
                value = (abs(T1.m - T2.m) < 10*eps(T1.m)+10*eps(T2.m)); 
            end
            value = all(value(:));
        end
        function value = ne(T1,T2) % overload ~=
            value = ~eq(T1,T2);
        end
        function T1 = plus(T1,T2) % overloaded + operator
            T1.m = T1.m + T2.m; 
        end
        function T1 = minus(T1,T2) % overloaded - operator
            T1.m = T1.m - T2.m; 
        end 
        function T1 = uplus(T1) % overloaded + operator (unitary)
            % pass
        end
        function T1 = uminus(T1) % overloaded - operator (unitary)
            T1.m = -T1.m; 
        end 
        function T1 = times(T1,T2) % overloaded .* (multiplication by scalar)
            if isa(T1,'anakin.tensor') % then b is not ten
                T1.m = T1.m .* T2; 
            else % a is not vector
                T2.m = T1 .* T2.m; 
                T1 = T2;
            end                
        end
        function T1 = mtimes(T1,T2) % overloaded * (multiplication by scalar or matrix)
            if isa(T1,'anakin.tensor') % then b is not vector
                T1.m = T1.m * T2; 
            else % a is not vector
                T2.m = T1 * T2.m; 
                T1 = T2;
            end
        end 
        function T1 = rdivide(T1,x) % overloaded ./ (division by scalar)
            T1.m = T1.m ./ x; 
        end
        function T1 = mrdivide(T1,x) % overloaded / (division by scalar or matrix)
            T1.m = T1.m / x; 
        end
        function T1 = ldivide(x,T1) % overloaded .\ (division by scalar)
            T1.m = x .\ T1.m; 
        end
        function T1 = mldivide(x,T1) % overloaded \ (division by scalar or matrix)
            T1.m = x \ T1.m; 
        end         
        function [l1,l2,l3,v1,v2,v3] = eigs(T) % like Matlab eigs function
           [V,D] = eigs(T.m);
           l1 = D(1,1);
           l2 = D(2,2);
           l3 = D(3,3);           
           v1 = anakin.vector(V(:,1));
           v2 = anakin.vector(V(:,2));
           v3 = anakin.vector(V(:,3));
        end
        function disp(T) % display
            disp('Tensor with canonical components:')
            disp(T.m)
        end
    end 
    methods % functionality
        function matrix = matrix(T,B1) % matrix of the tensor in B1
            if ~exist('B1','var')
                B1 = anakin.basis; % default frame
            end 
            matrix = B1.m' * T.m * B1.m;
            if isa(matrix,'sym')
                matrix = formula(simplify(matrix));
            end
        end 
    end
    methods % symbolic
        function dI = dt(T,B) % time derivative with respect to basis B. Requires sym vector that utlimately depends on a single variable t
            if ~exist('B','var')
                B = anakin.basis; % canonical vector basis
            end
            dI = anakin.tensor(diff(sym(T.matrix(B)),1),B);
        end
        function T_ = subs(T,variables,values) % particularize symbolic vector
            T_ = T;
            T_.c = double(subs(T.m,variables,values));
        end
    end
    methods % logical tests
        function isunitary = isunitary(T) % tensor is Hermitian (A' * A = eye(3))
            if isa(a.c,'sym') % symbolic inputs
                isunitary = isAlways(T.m' * T.m == eye(3),'Unknown','false'); % In case of doubt, false
            else % numeric input            
                isunitary = (abs(T.m' * T.m - eye(3))<eps(max(abs(T.m)))); 
            end 
            isunitary = all(isunitary);
        end 
    end
end








