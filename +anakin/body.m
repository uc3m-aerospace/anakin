%{
DESCRIPTION:
particle: class to model a rigid body. 

SYNTAX:
b0 = anakin.body();  % returns default object  
b  = anakin.body(<b>|(<<<mass>,inertia>,S|(<A>,<B>)>),<S1>);
where:
- <> denotes optional arguments
- | denotes alternative arguments
- b0 is the default bidy (mass 1, inertia eye(3), G located at the origin)
- b  is a particle  
- mass is a scalar
- inertia is a 2nd order tensor with the tensor of inertia at the center of
  mass
- S is a frame
- A is a point (center of mass) 
- S1 is a frame. If given, all previous input as relative to that frame

PROPERTIES:
* mass: the mass of the particle
* point: the point where the center of mass of the body is
* forces: a cell array with all the vector forces acting on the particle
* torques: a cell array with all the vector torques at G acting on the
  particle (must add manually the torques due to applied forces)

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
classdef body < anakin.frame & anakin.particle
    properties (Hidden = true, Access = protected) 
        I anakin.tensor = anakin.tensor([1,0,0;0,1,0;0,0,1]); % tensor of inertia of the object at the center of mass (canonical components)
    end
    properties 
        torques cell = {}; % cell array with all torques (vectors) acting on the object
    end
    methods % creation b  = anakin.body(<b|<<mass>,inertia>,S|(<A>,<B>)>,<S1>);
        function b = body(varargin) % constructor
            switch nargin
                case 0 % no arguments       
                    return;
                case 1 % body
                    b.M = varargin{1}.M;
                    b.I = varargin{1}.I;
                    b.v = varargin{1}.v;
                    b.b = varargin{1}.b;
                case 4 % mass, inertia, vector, basis
                    b.M = anakin.tensor(varargin{1});
                    b.I = anakin.tensor(varargin{2});
                    b.v = anakin.tensor(varargin{3});
                    b.b = anakin.basis(varargin{4});                    
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in body');
            end     
        end  
        function b = set.I(b,value) % on setting i
            b.I = anakin.tensor(value); 
        end 
        function b = set.torques(b,value) % on setting torques
            if ~iscell(value) % ensure cell array
                value = {value};
            end
            for j=1:length(value) % validate input
                if ~isa(value{j},'anakin.tensor') 
                    error('The torques must be supplied in a cell array of anakin.tensor vectors');
                end
            end
            b.torques = value; 
        end
    end 
    methods (Hidden = true) % overloads
        function value = eq(b1,b2) % overload ==
            value = (b1.M == b2.M) && (b1.I == b2.I) && (b1.v == b2.v) && (b1.b == b2.b);
        end
        function value = ne(b1,b2) % overload =~
            value = ~eq(b1,b2);
        end
        function disp(P) % display
            disp('Rigid body with mass:')
            disp(P.M.components)      
            disp('Inertia tensor at the center of mass:')
            disp(P.I.components)     
            disp('Coordinates of the center of mass:')
            disp(P.v.components)            
            disp('And basis with rotation matrix:')
            disp(P.b.matrix)                
        end
    end
    methods % general functionality  
        function frame = frame(b)
            frame = anakin.frame(b.v,b.b);
        end
        function inertia = inertia(b,O) % inertia tensor of the body with respect to point O
            if exist('O','var')
                r = b.center.pos - O.pos;
                inertia = b.I + b.M * (norm(r)^2 * anakin.tensor(eye(3)) - product(r,r));
            else
                inertia = b.I; % inertia of the body about its center of mass
            end            
        end
        function pointpos = pointpos(b,P) % position of a bodypoint
            pointpos = P.pos-b.pos;
        end
        function pointvel = pointvel(b,P,S1) % velocity of a bodypoint
            pointvel = b.pointpos(P).vel(S1);
        end
        function pointaccel = pointaccel(b,P,S1) % acceleration of a bodypoint
            pointaccel = b.pointpos(P).accel(S1);
        end
        function H = H(b,O,S1) % angular momentum about O in S1
            if ~exist('O','var')
                O = anakin.point; % default point
            end
            if exist('S1','var')
                H = cross(b.pos-O.pos, b.p(S1)) + b.I * b.omega(S1);
            else
                H = cross(b.pos-O.pos, b.p) + b.I * b.omega;
            end            
        end
        function T = T(b,S1) % kinetic energy in S1
            if exist('S1','var')
                vel = b.vel(S1).components;
                omega = b.omega(S1).components;
            else
                vel = b.vel.components;
                omega = b.omega.components;
            end
            T = b.M * norm(vel)^2 / 2 + omega * b.I * omega / 2; 
        end
        function eqs = equations(b,B1) % returns vector of equations of motion projected in basis B1
            p = b.p.dt;
            H = b.H.dt;
            F = anakin.tensor([0;0;0]); % allocate;
            M = anakin.tensor([0;0;0]); % allocate;
            for i=1:length(b.forces)
                F = F + b.forces{i};
                M = M + b.torques{i};
            end
            if ~exist('B1','var')
                B1 = anakin.basis;
            elseif isa(B1,'anakin.frame') 
                B1 = B1.basis; % extract basis
            end
            eqs = sym([0;0;0]);
            for i=1:3
                p_ = p * B1.e(i);
                F_ = F * B1.e(i);
                eqs(i) = (p_.components == F_.components);
            end                
            for i=1:3
                H_ = H * B1.e(i);
                M_ = M * B1.e(i);
                eqs(3+i) = (H_.components == M_.components);
            end      
        end
        function b_ = subs(b,variables,values) % particularize symbolic vector
            b_ = b;
            b_.M = b.M.subs(variables,values);
            b_.I = b.I.subs(variables,values);
            b_.v = b.v.subs(variables,values);
            b_.b = b.b.subs(variables,values);
        end
    end      
end




