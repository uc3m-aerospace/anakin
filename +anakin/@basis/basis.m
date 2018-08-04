%{
basis: class to define vector bases.

A basis can be created by passing the transformation matrix or the
vectors of the basis. By default, the code assumes the canonical vector
basis B0 is being used in the definition. Optionally, a different basis
for the definition can be given as a last argument. Called without
arguments, basis returns the canonical vector basis B0.
  
Common operators have been overloaded so that bases multiplied or
compared. 

METHODS:
* matrix: transformation matrix to another basis 
* i,j,k: returns the vectors of the basis (with respect to a chosen
  basis) 
* rotatex,rotatey,rotatez: create simply rotated frame about one
  coordinated axis of another basis
* omega: returns the angular velocity vector omega with respect to
  another basis (symbolic variables must be used)
* alpha: returns the angular acceleration vector alpha with respect to
  another basis (symbolic variables must be used)
* subs: takes values of the symbolic unknowns and returns a basis with
  purely numeric matrix (symbolic variables must be used)    
* isorthonormal, isrighthanded: checks the corresponding property and
  returns true or false  
* isproper: isorthonormal && isrighthanded
* plot: plots the basis with quiver, at a chosen position

MMM20180802
%}
classdef basis
    properties (Hidden = true, Access = protected)        
        m = [1,0,0;0,1,0;0,0,1]; % transformation matrix to the canonical vector basis   
    end 
    methods % creation
        function B = basis(varargin) % creator
            switch nargin
                case 0 % no arguments 
                    return;
                case 1 
                    if isa(varargin{1},'anakin.basis') % basis
                        B.m = varargin{1}.m;
                    else % matrix
                        B.m = varargin{1};
                    end
                case 2 
                    if isa(varargin{1},'anakin.basis') % relative basis, basis
                        B.m = varargin{2}.m * varargin{1}.m;
                    else % relative matrix, basis
                        B.m = varargin{2}.m * varargin{1};
                    end
                case 3 % i,j,k (vectors or columns)
                    i = anakin.vector(varargin{1}); % ensure vector
                    j = anakin.vector(varargin{2});
                    k = anakin.vector(varargin{3});
                    B.m = [i.components,j.components,k.components];
                case 4 % relative i,j,k (vectors or columns) and basis
                    i = anakin.vector(varargin{1},varargin{4});
                    j = anakin.vector(varargin{2},varargin{4});
                    k = anakin.vector(varargin{3},varargin{4});
                    B.m = [i.components,j.components,k.components];
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in basis');
            end     
        end 
        function B = set.m(B,value) % on setting m
            B.m = reshape(value,3,3);
            if isa(B.m,'sym') % symbolic input
                B.m = formula(simplify(B.m)); % simplify and force sym rather than symfun to allow indexing
            end
        end     
    end
    methods % overloads 
        function B3 = mtimes(B1,B2) % overloaded * (multiplication of two bases)
            B3 = B1;
            B3.m = B1.m * B2.m;
        end 
        function B3 = mrdivide(B1,B2) % overloaded /
            B3 = B1;
            B3.m = B1.m / B2.m;
        end 
        function B3 = mldivide(B1,B2) % overloaded \
            B3 = B1;
            B3.m = B1.m \ B2.m;
        end 
        function value = eq(B1,B2) % overload ==
            if isa(B1.m,'sym') || isa(B1.m,'sym') % symbolic inputs
                value = isAlways(B1.m==B2.m,'Unknown','false'); % In case of doubt, false
            else % numeric input
                value = (abs(B1.m - B2.m)<eps(B1.m)+eps(B2.m)); 
            end
            value = all(value(:));
        end
        function value = ne(B1,B2) % overload ~=
            value = ~eq(B1,B2);
        end
    end
    methods % functionality
        function matrix = matrix(B,B1) % transformation matrix to another basis
            if ~exist('B1','var')
                matrix = B.m; % if no basis is given, use the canonical vector basis
            else
                matrix = B1.m' * B.m;
            end
            if isa(matrix,'sym')
                matrix = simplify(matrix);
            end
        end
        function i = i(B) % vector i of the basis
             i = anakin.vector([1;0;0],B);
        end
        function j = j(B) % vector j of the basis
             j = anakin.vector([0;1;0],B);
        end
        function k = k(B) % vector k of the basis
             k = anakin.vector([0;0;1],B);
        end 
        function Bx = rotatex(B,angle) % returns rotated basis about x axis of B by angle
            Bx = B;
            Bx.m = B.m * [1,0,0;0,cos(angle),-sin(angle);0,sin(angle),cos(angle)];            
        end
        function By = rotatey(B,angle) % returns rotated basis about y axis of B by angle
            By = B;
            By.m = B.m * [cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
        end
        function Bz = rotatez(B,angle) % returns rotated basis about z axis of B by angle
            Bz = B;
            Bz.m = B.m * [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
        end        
        function omega = omega(B,B1) % Returns the symbolic angular velocity vector with respect to B1
            omega = anakin.vector([dot(B.k,B.j.dt); dot(B.i,B.k.dt); dot(B.j,B.i.dt)],B); % If B1 is not given, assume the canonical vector basis B0
            if exist('B1','var') % If B1 is given, correct previous value
                omega = omega - B1.omega; 
            end 
        end
        function alpha = alpha(B,B1) % Returns the symbolic angular acceleration vector with respect to B1
            alpha = B.omega.dt; % If B1 is not given, assume the canonical vector basis B0
            if exist('B1','var') % If B1 is given, correct previous value
                alpha = alpha - cross(B1.omega,B.omega) - B1.alpha; 
            end
        end
        function B_ = subs(B,variables,values) % particularize symbolic basis
            B_ = B;
            B_.m = double(subs(B.m,variables,values));
        end
    end
    methods % logical tests
        function isorthonormal = isorthonormal(B) % all vectors are unitary and mutually orthogonal
            isorthonormal = (B.m' * B.m == eye(3));
            if isa(isorthonormal,'sym')
                isorthonormal = isAlways(isorthonormal,'Unknown','false'); % In case of doubt, false
            end
            isorthonormal = all(isorthonormal(:));
        end    
        function isrighthanded = isrighthanded(B) % basis is righthanded
            isrighthanded = (det(B.m) > 0);
            if isa(isrighthanded,'sym')
                isrighthanded = isAlways(isrighthanded,'Unknown','false'); % In case of doubt, false
            end
        end  
        function isproper = isproper(B) % orthonormal and righthanded
            isproper = B.isorthonormal && B.isrighthanded;
        end    
    end
    methods % plotting
        function h = plot(B,varargin) % plot. First argument in varargin must be the O vector, if any
            if mod(nargin,2) == 1 % no origin vector is given
                O = anakin.vector; % null vector
            else
                O = varargin{1};
                varargin = varargin(2:end);
            end
            cc = O.components;
            mm = B.m;
            hold on            
            h = quiver3([cc(1),cc(1),cc(1)],[cc(2),cc(2),cc(2)],[cc(3),cc(3),cc(3)],...
            mm(1,:),mm(2,:),mm(3,:),0,'color','k');
            hold off
            if ~isempty(varargin)
                set(h,varargin{:}); % set options stored in varargin
            end
        end
    end
end
