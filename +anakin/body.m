%{
DESCRIPTION:
body: class to model a rigid body. 

SYNTAX:
b0 = anakin.body();  % returns default object  
b  = anakin.body(<b>,<S1>);
b  = anakin.body(<mass>,<inertia>,<S>,<S1>);
b  = anakin.body(<mass>,<inertia>,<A|a|c>,<B>,<S1>);
where:
- <> denotes optional arguments
- | denotes alternative arguments
- b0 is the default bidy (mass 1, inertia eye(3), G located at the origin)
- b  is a particle  
- S is a frame
- mass is a scalar or number
- inertia is a 2nd order tensor with the tensor of inertia at G in basis B
- A is a point (center of mass G)
- c is an array with the Cartesian coordinates of the origin
- B  is a basis
- S1 is a frame. If given, all previous input as relative to that frame
 
METHODS: 
* p: linear momentum in a given reference frame
* H: angular momentum about a point in a given reference frame
* T: kinetic energy in a given reference frame
* equations: returns a vector of (symbolic) equations, m*a = F, dH/dt = M,
  projected along the vectors of one basis
* inertia: tensor of inertia about a point
* subs: takes values of the symbolic unknowns and returns a particle
  object which is purely numeric

To add: instantaneous axis of rotation and slip, slip velocity, velocity of
a point, acceleration of a point, add inertia tensor given at a point,
facilities to add forces at a point, torques, visualize them, frame
AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef body < anakin.frame
    properties (Hidden = true, Access = protected) 
        M anakin.tensor = anakin.tensor(1); % mass of the object
        I0 anakin.tensor = anakin.tensor([0,0,0;0,0,0;0,0,0]); % tensor of inertia of the body at the center of mass in the canonical vector basis
    end 
    methods % creation b  = anakin.body(<b|<<mass>,inertia>,S|(<A>,<B>)>,<S1>);
        function b = body(varargin) % constructor
            if isempty(varargin) % Default body
                return;
            elseif length(varargin) > 1 && isa(varargin{end},'anakin.frame') % Last argin is frame
                S1 = varargin{end};
                varargin = varargin(1:end-1); 
            else % No frame is provided; use default
                S1 = anakin.frame;
            end 
            b.v = S1.v;
            b.m = S1.m;
            for i = 1:length(varargin) % later inputs overwrite former inputs
                temp = varargin{i};
                if isa(temp,'anakin.body')
                    b.M = temp.M;
                    b.I0 = anakin.tensor(temp.I0,S1);
                    b.v = anakin.tensor(S1.v.components + S1.m * temp.v.components);
                    b.m = temp.matrix(S1);
                elseif isa(temp,'anakin.particle')
                    b.M = temp.M;
                    b.v = anakin.tensor(S1.v.components + S1.m * temp.v.components);
                elseif isa(temp,'anakin.frame')
                    b.v = anakin.tensor(S1.v.components + S1.m * temp.v.components);
                    b.m = temp.matrix(S1);
                elseif isa(temp,'anakin.point')
                    b.v = anakin.tensor(S1.v.components + S1.m * temp.v.components); 
                elseif isa(temp,'anakin.basis')
                    b.m = temp.matrix(S1);
                elseif isa(temp,'anakin.tensor')
                    if temp.ndims == 0 % assume it is M
                        b.M = temp;
                    elseif temp.ndims == 1 % assume it is v
                        b.v = anakin.tensor(S1.v.components + S1.m * temp.components);
                    elseif temp.ndims == 2 % assume it is I0
                        b.I0 = anakin.tensor(temp.components,S1);
                    else
                        error('Cannot take tensors of order higher than 2 as inputs');
                    end
                else % Array
                    if numel(temp) == 1 % assume it is M
                        b.M = temp;
                    elseif length(temp(:,1)) > 1 && length(temp(1,:)) > 1 % assume it is I0. No B matrix is allowed this way.
                        b.I0 = anakin.tensor(temp,S1);
                    else % assume it is v
                        v_ = anakin.tensor(temp);
                        b.v = anakin.tensor(S1.v.components + S1.m * v_.components);
                    end
                end
            end 
        end
        function b = set.M(b,value) % on setting M
            b.M = anakin.tensor(value); 
        end 
        function b = set.I0(b,value) % on setting i
            b.I0 = anakin.tensor(value); 
        end  
    end 
    methods (Hidden = true) % overloads
        function value = eq(b1,b2) % overload ==  
            I01 = anakin.tensor(b1.I0,b1.basis);
            I02 = anakin.tensor(b2.I0,b2.basis); 
            value = (b1.M == b2.M) && (I01 == I02) && (b1.v == b2.v);
        end
        function value = ne(b1,b2) % overload =~
            value = ~eq(b1,b2);
        end
        function disp(b) % display
            disp('Rigid body with mass:')
            disp(b.M.components)      
            disp('Inertia tensor at the center of mass in body basis:')
            I = b.I0.components(b);
            try
                I = double(I);
            catch
                % pass
            end
            disp(I)
            disp('Coordinates of the center of mass:')
            disp(b.v.components)            
            disp('And basis with rotation matrix:')
            disp(b.m)                
        end
    end
    methods % general functionality   
        function mass = mass(b) % mass of the particle
            mass = b.M;
        end
        function p = p(b,S1) % linear momentum in S1
            if exist('S1','var')
                p = b.M*b.vel(S1);
            else
                p = b.M*b.vel;
            end            
        end
        function H = H(b,O,S1) % angular momentum about O in S1
            if ~exist('O','var')
                O = anakin.point; % default point
            end
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end  
            H = cross(b.pos-O.pos, b.p(S1)) + b.I0 * b.omega(S1); 
        end
        function T = T(b,S1) % kinetic energy in S1
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end  
            vel = b.vel(S1);
            omega = b.omega(S1);
            T = b.M * norm(vel)^2 / 2 + omega * b.I0 * omega / 2; 
        end 
        function I = I(b,O) % inertia tensor of the body with respect to point O in canonical vector basis
            if ~exist('O','var')
                O = anakin.point; % default point is the origin
            end
            r = b.v - O.v;
            I = b.I0 + b.M * (norm(r)^2 * eye(3) - product(r,r)); 
        end
        function b_ = subs(b,variables,values) % particularize symbolic vector
            b_ = b;
            b_.M = b.M.subs(variables,values);
            b_.I0 = b.I0.subs(variables,values);
            b_.v = b.v.subs(variables,values); 
            b_.m = subs(b.m,variables,values);
            try
                b_.m = double(b_.m);
            catch
                % pass
            end
        end
    end      
end




