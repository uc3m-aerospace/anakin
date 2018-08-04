%{
vector: class to define a 3-vector.

A vector can be created by passing coordinates as a column array or as
three independent values. By default, the code assumes the canonical
vector basis B0 is being used in the definition. Optionally, a different
basis for the definition can be given as a last argument. Called without
arguments, the null vector is returned. 

Common operators have been overloaded so that vectors can be added,
subtracted, dotted, crossed, compared, etc.

METHODS:
* components: returns the components of the vector in a chosen basis
* x,y,z: returns individual components in a chosen basis
* isunitary, isperpendicular, isparallel: checks for the corresponding
  property and returns true or false    
* magnitude, dir: returns the magnitude and unit vector of a vector
* dt: returns the time derivative of a vector wrt a chosen basis
  (symbolic variables must be used)  
* subs: takes values of the symbolic unknowns and returns a vector with
  purely numeric coordinates (symbolic variables must be used)   
* plot: plots the vector with quiver, at a chosen position

MMM20180802
%}
classdef vector
    properties (Hidden = true, Access = protected) 
        c = [0;0;0]; % components of the vector in the canonical vector basis
    end 
    methods           
        function a = vector(varargin) % creator 
            switch nargin
                case 0 % no arguments: use defaults
                    return;
                case 1 % argument is: components in column vector in the canonical vector basis
                    if isa(varargin{1},'anakin.vector') % catch case where vector is called with a vector already
                        a = varargin{1};
                    else % input is a column array
                        a.c = varargin{1};
                    end
                case 2 % arguments are: components in column vector, basis B
                    a.c = varargin{2}.matrix * reshape(varargin{1},3,1);
                case 3 % arguments are: x, y, z
                    a.c = [varargin{1};varargin{2};varargin{3}];
                case 4 % arguments are: x, y, z, basis B 
                    a.c = varargin{4}.matrix * [varargin{1};varargin{2};varargin{3}];
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in vector');
            end       
        end
        function a = set.c(a,value)
            a.c = reshape(value,3,1); % Force column 
            if isa(a.c,'sym') % symbolic input
                a.c = formula(simplify(a.c)); % simplify and force sym rather than symfun to allow indexing into c
            end
        end
    end
    methods 
        function a = plus(a,b) % overloaded + operator 
            a.c = a.c + b.c; 
        end
        function a = minus(a,b) % overloaded - operator 
            a.c = a.c - b.c; 
        end 
        function a = uplus(a) % overloaded + operator (unitary)
            % pass
        end
        function a = uminus(a) % overloaded - operator (unitary)
            a.c = -a.c; 
        end 
        function a = times(a,b) % overloaded .* (multiplication by scalar)
            if isa(a,'anakin.vector') % then b is scalar
                a.c = a.c.*b; 
            else 
                b.c = a.*b.c; 
                a = b;
            end                
        end
        function a = mtimes(a,b) % overloaded * (multiplication by scalar or matrix)
            if isa(a,'anakin.vector') % then b is scalar or matrix
                a.c = a.c*b; 
            else 
                b.c = a*b.c; 
                a = b;
            end
        end 
        function a = rdivide(a,x) % overloaded ./ (division by scalar)
            a.c = a.c./x; 
        end
        function a = mrdivide(a,x) % overloaded / (division by scalar or matrix)
            a.c = a.c/x; 
        end
        function a = ldivide(x,a) % overloaded .\ (division by scalar)
            a.c = x.\a.c; 
        end
        function a = mldivide(x,a) % overloaded \ (division by scalar or matrix)
            a.c = x\a.c; 
        end
        function value = eq(a,b) % overload ==
            if isa(a.c,'sym') || isa(b.c,'sym') % symbolic inputs
                value = isAlways(a.c==b.c,'Unknown','false'); % In case of doubt, false
            else % numeric input            
                value = (abs(a.c - b.c)<eps(a.c)+eps(b.c)); 
            end
            value = all(value(:));
        end
        function value = ne(a,b) % overload ~=
            value = ~eq(a,b);
        end
        function value = dot(a,b) % dot product of two real vectors
            value = dot(a.c,b.c);
            if isa(value,'sym')
                value = formula(simplify(value)); % simplify and force sym rather than symfun to allow indexing
            end
        end
        function value = norm(a) % 2 norm of a real vector
            value = dot(a,a);
            if isa(value,'sym')
                value = formula(simplify(value)); % simplify and force sym rather than symfun to allow indexing
            end
        end
        function value = cross(a,b) % cross product
            value = anakin.vector(cross(a.c,b.c)); 
        end        
    end
    methods
        function dir = dir(a) % returns unit vector along a
            dir = a/norm(a);
        end
        function magnitude = magnitude(a) % returns magnitude of a (alias for norm)
            magnitude = norm(a);
        end
    end
    methods        
        function components = components(a,B) % return column of components of a in basis B
            if ~exist('B','var')
                components = a.c; % if no basis is given, use the canonical vector basis
            else
                B0 = anakin.basis; % canonical vector basis
                components = B0.matrix(B) * a.c;
            end
            if isa(components,'sym')
                components = simplify(components);
            end
        end
        function x = x(a,B) % returns a single component in basis B
            if ~exist('B','var')
                B = anakin.basis; % canonical vector basis
            end
            components = a.components(B);
            x = components(1);
        end
        function y = y(a,B)
            if ~exist('B','var')
                B = anakin.basis; % canonical vector basis
            end
            components = a.components(B);
            y = components(2);
        end
        function z = z(a,B)
            if ~exist('B','var')
                B = anakin.basis; % canonical vector basis
            end
            components = a.components(B);
            z = components(3);
        end 
        function a_ = subs(a,variables,values) % particularize symbolic vector
            a_ = anakin.vector(double(subs(a.c,variables,values)));
        end
        function da = dt(a,B) % time derivative with respect to basis B. Requires sym vector that utlimately depends on a single variable t
            if ~exist('B','var')
                B = anakin.basis; % canonical vector basis
            end
            da = anakin.vector(diff(sym(a.components(B)),1),B);
        end
    end
    methods
        function isunitary = isunitary(a) % vector is unitary
            isunitary = (dot(a,a)==1);
            if isa(isunitary,'sym')
                isunitary = isAlways(isunitary,'Unknown','false'); % In case of doubt, false
            end
        end
        function isperpendicular = isperpendicular(a,b) % the two vectors are perpendicular
            isperpendicular = (dot(a,b)==0);
            if isa(isperpendicular,'sym')
                isperpendicular = isAlways(isperpendicular,'Unknown','false'); % In case of doubt, false
            end
        end
        function isparallel = isparallel(a,b) % the two vectors are parellel
            isparallel = (cross(a,b)==anakin.vector(0,0,0));
            if isa(isparallel,'sym')
                isparallel = isAlways(isparallel,'Unknown','false'); % In case of doubt, false
            end
        end     
    end
    methods
        function h = plot(v,varargin) % plot. First argument in varargin must be the O vector, if any
            if mod(nargin,2) == 1 % no origin vector is given
                O = anakin.vector; % null vector
            else
                O = varargin{1};
                varargin = varargin(2:end);
            end 
            hold on
            h = quiver3(O.c(1),O.c(2),O.c(3),v.c(1),v.c(2),v.c(3),0,'color','k');
            hold off
            if ~isempty(varargin)
                set(h,varargin{:}); % set options stored in varargin
            end
        end
    end
end








