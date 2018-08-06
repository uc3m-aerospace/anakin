%{
basis: class to define orthonormal, right-handed vector bases.

The class constructor accepts the following call types:
- B0 = anakin.basis(); % where: B0 is the canonical vector basis
- B = anakin.basis(B); % (convert to basis class)
- B = anakin.basis(m); % where: m is rotation matrix from B0 to B
- B = anakin.basis(m,B1); % where: m is relative rotation matrix from B1
                          % to B 
- B = anakin.basis(i,j,k); % where: i,j,k are vector objects
- B = anakin.basis(ic,jc,kc); % where: ic,jc,kc are column components in
                              % B0
- B = anakin.basis(ric,rjc,rkc,B1); % where: ric,rjc, rkc are relative
                                    % column components, and B1 is a
                                    % vector basis  
  
METHODS:
* matrix: transformation matrix to another basis 
* i,j,k: returns the vectors of the basis (with respect to a chosen
  basis) 
* axis, angle: returns unit vector and angle of rotation of B wrt B1
* quaternions: returns the quaternions with respect to another basis
* rotatex,rotatey,rotatez: create simply rotated frame about one
  coordinated axis of another basis
* omega, alpha: returns the angular velocity/acceleration vector omega
  with respect to another basis (symbolic variables must be used)
* subs: takes values of the symbolic unknowns and returns a basis with
  purely numeric matrix (symbolic variables must be used)    
* isorthonormal, isrighthanded: checks the corresponding property and
  returns true or false  
* plot: plots the basis with quiver, at a chosen position

MMM20180802
%}
classdef basis
    properties (Hidden = true, Access = protected)        
        m = [1,0,0;0,1,0;0,0,1]; % transformation matrix: [a(in B0)] = m * [a(in B)]. Or equivalently: the rotation matrix to go from B0 to B
    end 
    methods % creation
        function B = basis(varargin) % constructor
            switch nargin
                case 0 % no arguments 
                    return; 
                case 1
                    if isa(varargin{1},'anakin.basis') % basis
                        B.m = varargin{1}.m;
                    else % matrix
                        B.m = varargin{1};
                    end
                case 2 % relative matrix, basis
                    B.m = varargin{2}.m * varargin{1};                     
                case 3 % i,j,k (vectors or component columns)
                    B.m = [anakin.vector(varargin{1}).components,...
                           anakin.vector(varargin{2}).components,...
                           anakin.vector(varargin{3}).components];
                case 4 % relative i,j,k (component columns) and basis
                    B.m = [anakin.vector(varargin{1},varargin{4}).components,...
                           anakin.vector(varargin{2},varargin{4}).components,...
                           anakin.vector(varargin{3},varargin{4}).components];
                otherwise
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
        function value = eq(B1,B2) % overload ==
            if isa(B1.m,'sym') || isa(B1.m,'sym') % symbolic inputs
                value = isAlways(B1.m==B2.m,'Unknown','false'); % In case of doubt, false
            else % numeric input
                value = (abs(B1.m - B2.m)<eps(max(abs(B1.m(:))))+eps(max(abs(B2.m(:))))); 
            end
            value = all(value(:));
        end
        function value = ne(B1,B2) % overload ~=
            value = ~eq(B1,B2);
        end
    end
    methods % functionality
        function matrix = matrix(B,B1) % transformation matrix to another basis: [a(in B1)] = m * [a(in B)]
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
        function axis = axis(B,B1) % rotation axis unit vector from B1
            if ~exist('B1','var')
                B1 = anakin.basis; % if no basis is given, use the canonical vector basis
            end
            mm = B.matrix(B1);
            axis = anakin.vector([mm(2,3)-mm(3,2);mm(3,1)-mm(1,3);mm(1,2)-mm(2,1)],B1).dir; % fails if rotation angle is 0 or 180 deg
        end 
        function angle = angle(B,B1) % angle of rotation from B1
            if ~exist('B1','var')
                B1 = anakin.basis; % if no basis is given, use the canonical vector basis
            end
            mm = B.matrix(B1);
            angle = acos((trace(mm)-1)/2);
            if isa(angle,'sym') % symbolic input
                angle = formula(simplify(angle)); % simplify and force sym rather than symfun
            end

        end
        function quaternions = quaternions(B,B1) % quaternions of rotation from B1. Fails when rotation angle is 180 deg
            if ~exist('B1','var')
                B1 = anakin.basis; % if no basis is given, use the canonical vector basis
            end
            mm = B.matrix(B1);
            quaternions(4) = sqrt(trace(mm)+1)/2; % scalar term q4
            quaternions(1) = -(mm(2,3)-mm(3,2))/(4*quaternions(4)); % q1
            quaternions(2) = -(mm(3,1)-mm(1,3))/(4*quaternions(4)); % q2
            quaternions(3) = -(mm(1,2)-mm(2,1))/(4*quaternions(4)); % q3
            quaternions = reshape(quaternions,4,1); % Force column
            if isa(quaternions,'sym') % symbolic input
                quaternions = formula(simplify(quaternions)); % simplify and force sym rather than symfun to allow indexing
            end
        end
        function euler = euler(B,type,B1) % Euler angles of rotation from B1. Fails depending on the value of the intermediate angle: symmetric Euler angles fail for theta2 = 0,180 deg. Asymmetric Euler angles fail for theta2 = 90,270 deg 
            % type: a text string like '313' or '123' indicating the
            % intrinsic axes of rotation
            if ~exist('B1','var')
                B1 = anakin.basis; % if no basis is given, use the canonical vector basis
            end
            tbi(B,type,B1);
            euler = 'to be implemented';
            error(euler);
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
