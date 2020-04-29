%{
DESCRIPTION:
basis: class to define orthonormal, right-handed vector bases.

SYNTAX:
B = anakin.basis();  % returns default object (canonical vector basis)
B = anakin.basis(B,<B1>);
B = anakin.basis(m,<B1>);
B = anakin.basis(v,v,v,<B1>);
B = anakin.basis(c,c,c,<B1>);
B = anakin.basis(q,<B1>);
B = anakin.basis(axis,angle,<B1>);
where:
- <> denotes optional arguments 
- B  is a basis
- m  is a square rotation matrix (square array)
- v is a vector (1st-order tensor)
- c is an array with the vector components 
- q are quaternions (with the scalar last). Only used for 3D bases
- axis is the unit vector of the axis of rotation
- angle is angle of rotation about axis
- B1 is a basis. If given, all previous inputs are relative to that basis

METHODS:
* spacedim: returns dimensionality of space
* matrix: transformation matrix to another basis 
* e: returns the i-th vector of the basis
* isisorthonormal: true if the transformation matrix is orthonormal
* subs: takes values of the symbolic unknowns and returns a basis with
  purely numeric matrix (symbolic variables must be used)    
* plot: plots vectors of basis

(three dimensional space only):
* rotaxis, rotangle: returns unit vector and angle of rotation of B wrt
  B1 
* quaternions: returns the quaternions with respect to another basis
* euler: returns the euler angles (of a chosen type) with respect to
  another basis
* rotatex,rotatey,rotatez: create simply rotated frame about one
  coordinated axis of another basis
* isrighthanded: true if the basis is right handed
* omega, alpha: returns the angular velocity/acceleration vector omega
  with respect to another basis (symbolic variables must be used)

AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef basis
    properties (Hidden = true)
        m = [1,0,0;0,1,0;0,0,1]; % transformation matrix: [a(in B0)] = m * [a(in B)]. Or equivalently: the rotation matrix to go from B to B0
    end 
    methods % creation
        function B = basis(varargin) % constructor
            for i = 1:length(varargin)
               if isa(varargin{i},'sym')
                   varargin{i} = formula(varargin{i}); % enforce formula to allow indexing
               end
               if isa(varargin{i},'anakin.frame')
                   varargin{i} = varargin{i}.basis; % take only basis
               end
            end  
            switch nargin
                case 0 % no arguments 
                    return;    
                case 1
                    B.m = anakin.basis(varargin{1},anakin.basis).m;              
                case 2 
                    if isa(varargin{end},'anakin.basis') 
                        if isa(varargin{1},'anakin.basis') % relative basis, basis
                            B.m = varargin{2}.m * varargin{1}.m; 
                        elseif length(varargin{1}) == 4 % relative quaternions, basis
                            qq = varargin{1}; % quaternions with the scalar component last
                            mm = [qq(4)^2+qq(1)^2-qq(2)^2-qq(3)^2,     2*(qq(1)*qq(2)-qq(4)*qq(3)),     2*(qq(1)*qq(3)+qq(4)*qq(2)); % matrix whose columns are the components of the ijk vectors of B expressed in B1
                                      2*(qq(1)*qq(2)+qq(4)*qq(3)), qq(4)^2-qq(1)^2+qq(2)^2-qq(3)^2,     2*(qq(2)*qq(3)-qq(4)*qq(1));
                                      2*(qq(1)*qq(3)-qq(4)*qq(2)),     2*(qq(2)*qq(3)+qq(4)*qq(1)), qq(4)^2-qq(1)^2-qq(2)^2+qq(3)^2];
                            B.m = varargin{2}.m * mm;
                        else % relative matrix, basis
                            if ismatrix(varargin{1}) && length(varargin{1}(:,1)) == length(varargin{1}(1,:))
                                B.m = varargin{2}.m * varargin{1};
                            else
                                error('The matrix is not square')
                            end
                        end
                    else 
                        B.m = anakin.basis(varargin{1},varargin{2},anakin.basis).m; 
                    end 
                case 3 
                    if isa(varargin{end},'anakin.basis') % relative axis, relative angle, basis
                        axis = anakin.tensor(varargin{1},varargin{3}).components;
                        angle = varargin{2};
                        c = cos(angle);
                        s = sin(angle);
                        C = 1-c;                        
                        B.m = [              axis(1)^2*C+c, axis(1)*axis(2)*C-axis(3)*s, axis(1)*axis(3)*C+axis(2)*s;
                               axis(1)*axis(2)*C+axis(3)*s,               axis(2)^2*C+c, axis(2)*axis(3)*C-axis(1)*s;
                               axis(1)*axis(3)*C-axis(2)*s, axis(2)*axis(3)*C+axis(1)*s,               axis(3)^2*C+c];
                    else 
                        B.m = anakin.basis(varargin{1},varargin{2},varargin{3},anakin.basis).m;  
                    end
                case 4 % relative i,j,k (component columns) and basis
                    B.m = [anakin.tensor(varargin{1},varargin{4}).components,...
                           anakin.tensor(varargin{2},varargin{4}).components,...
                           anakin.tensor(varargin{3},varargin{4}).components];                
                otherwise
                    error('Wrong number of arguments in basis');
            end  
        end 
        function B = set.m(B,value) % on setting m
            if isa(value,'sym') % symbolic input
                value = formula(simplify(value)); % simplify and force sym rather than symfun to allow indexing
            end
            B.m = value;
            try
                B.m = double(B.m);
            catch
                % pass
            end
        end     
    end
    methods (Hidden = true) % overloads 
        function value = eq(B1,B2) % overload ==
            if isa(B1.m,'sym') || isa(B1.m,'sym') % symbolic inputs
                value = isAlways(B1.m==B2.m,'Unknown','false'); % In case of doubt, false
            else % numeric input
                value = (abs(B1.m - B2.m) < 10*eps(B1.m)+10*eps(B2.m)); 
            end
            value = all(value(:));
        end
        function value = ne(B1,B2) % overload ~=
            value = ~eq(B1,B2);
        end 
        function disp(B) % display
            disp('Basis with rotation matrix:')
            disp(B.m)
        end
    end
    methods % general functionality
        function spacedim = spacedim(B) % number of dimensions of space
            spacedim = length(B.m(:,1));            
        end
        function matrix = matrix(B,B1) % transformation matrix to another basis: [a(in B1)] = m * [a(in B)]
            if ~exist('B1','var')
                matrix = B.m; % if no basis is given, use the canonical vector basis
            else
                matrix = B1.m \ B.m;
                if isa(matrix,'sym')
                    matrix = formula(simplify(matrix));
                end
            end            
        end        
        function e = e(B,i) % i-th vector of the basis
             e = anakin.tensor(B.m(:,i));
        end
        function isorthonormal = isorthonormal(B) % all vectors are unitary and mutually orthogonal
            if isa(B.m,'sym') % symbolic inputs
                isorthonormal = isAlways(B.m' * B.m == eye(B.spacedim),'Unknown','false'); % In case of doubt, false
            else % numeric input            
                isorthonormal = (abs(B.m' * B.m - eye(B.spacedim))<eps(max(abs(B.m(:))))); 
            end 
            isorthonormal = all(isorthonormal(:));
        end    
        function B_ = subs(B,variables,values) % particularize symbolic basis
            B_ = B;
            B_.m = subs(B.m,variables,values);
            try
                B_.m = double(B_.m);
            catch
                % pass
            end
        end
        function h = plot(B,varargin) % plot vectors of basis
            h = []; % initialize
            for i = 1:B.spacedim
                v = B.e(i);
                h = [h; v.plot(varargin{:})]; %#ok<AGROW>
            end
        end
    end
    methods % 3-dimensional-space functionality
        function axis = rotaxis(B,B1) % rotation axis unit vector from B1
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            if ~exist('B1','var')
                B1 = anakin.basis(eye(B.spacedim)); % if no basis is given, use the canonical vector basis
            end
            mm = B.matrix(B1);
            axis = anakin.tensor([mm(3,2)-mm(2,3);mm(1,3)-mm(3,1);mm(2,1)-mm(1,2)],B1).unitvector; % fails if rotation angle is 0 or 180 deg
        end 
        function angle = rotangle(B,B1) % angle of rotation from B1
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            if ~exist('B1','var')
                B1 = anakin.basis(eye(B.spacedim)); % if no basis is given, use the canonical vector basis
            end
            mm = B.matrix(B1);
            angle = acos((trace(mm)-1)/2);
            if isa(angle,'sym') % symbolic input
                angle = formula(simplify(angle)); % simplify and force sym rather than symfun
            end
        end
        function quaternions = quaternions(B,B1) % quaternions of rotation from B1. Fails when rotation angle is 180 deg
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            if ~exist('B1','var')
                B1 = anakin.basis(eye(B.spacedim)); % if no basis is given, use the canonical vector basis
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
        function euler = euler(B3,type,B1) % Euler angles of chosen type from B1. Fails depending on the value of the intermediate angle: symmetric Euler angles fail for theta2 = 0,180 deg. Asymmetric Euler angles fail for theta2 = 90,270 deg 
            if B3.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            if ~exist('B1','var')
                B1 = anakin.basis(eye(B3.spacedim)); % if no basis is given, use the canonical vector basis
            end
            if ~exist('type','var') % type: a vector like [3,1,3] or [1,2,3] indicating the intrinsic axes of rotation
                type = [3,1,3];            
            end
            m0 = eye(3); % this is B0 here!
            m3 = B3.matrix(B1);
            
            one = anakin.tensor(m0(:,type(1))); % first rotation direction
            three = anakin.tensor(m3(:,type(3))); % third rotation direction 
            if type(1)==type(3) % symmetric Euler angles
                euler(2) = one.angle(three).components;
                two = cross(one,three); 
            else % asymmetric Euler angles 
                even = det(m0(:,type)); % get even/odd of permutation of type 
                euler(2) = even * asin(dot(one.components,three.components));
                two = -even * cross(one,three);
            end 
            two = two.unitvector; % second rotation direction, normalized
            temp0 = anakin.tensor(m0(:,type(2)));
            euler(1) = temp0.angle(two,one).components;
            euler(3) = two.angle(anakin.tensor(m3(:,type(2))),three).components; 
            
            euler = reshape(euler,1,3); % Force row
            if isa(euler,'sym') % symbolic input
                euler = formula(simplify(euler)); % simplify and force sym rather than symfun to allow indexing
            end 
        end
        function Bx = rotatex(B,angle) % returns rotated basis about x axis of B by angle
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            Bx = anakin.basis;
            Bx.m = B.m * [1,0,0;0,cos(angle),-sin(angle);0,sin(angle),cos(angle)];            
        end
        function By = rotatey(B,angle) % returns rotated basis about y axis of B by angle
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            By = anakin.basis;
            By.m = B.m * [cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
        end
        function Bz = rotatez(B,angle) % returns rotated basis about z axis of B by angle
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            Bz = anakin.basis;
            Bz.m = B.m * [cos(angle),-sin(angle),0;sin(angle),cos(angle),0;0,0,1];
        end                
        function isrighthanded = isrighthanded(B) % basis is righthanded
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            isrighthanded = (det(B.m) > 0);
            if isa(isrighthanded,'sym')
                isrighthanded = isAlways(isrighthanded,'Unknown','false'); % In case of doubt, false
            end
        end    
        function omega = omega(B,B1) % Returns the symbolic angular velocity vector with respect to B1
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            omega = anakin.tensor([dot(B.m(:,3),B.e(2).dt.components); 
                                   dot(B.m(:,1),B.e(3).dt.components); 
                                   dot(B.m(:,2),B.e(1).dt.components)],B);
            if exist('B1','var') % If B1 is given, correct previous value
                omega = omega - B1.omega; 
            end 
        end
        function alpha = alpha(B,B1) % Returns the symbolic angular acceleration vector with respect to B1
            if B.spacedim ~= 3
                error('This functionality is only available for bases in 3D space');
            end
            alpha = B.omega.dt; % If B1 is not given, assume the canonical vector basis B0
            if exist('B1','var') % If B1 is given, correct previous value
                alpha = alpha - cross(B1.omega,B.omega) - B1.alpha; 
            end
        end         
    end
end
