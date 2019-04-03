%{
DESCRIPTION:
particle: class to model a point particle. 

SYNTAX:
P0 = anakin.particle();  % returns default object  
P  = anakin.particle(<P>|(<<mass>,A|a>),<S1>);
where:
- <> denotes optional arguments
- | denotes alternative arguments
- P0 is the default particle (mass 1, located at the origin)
- P  is a particle  
- mass is a scalar
- A is a point
- a is a vector 
- S1 is a frame. If given, all previous input as relative to that frame

PROPERTIES:
* mass: the mass of the particle
* point: the point where the particle is
* forces: a cell array with all the vector forces acting on the particle

METHODS: 
* mass: the mass of the particle
* p: linear momentum in a given reference frame
* H: angular momentum about a point in a given reference frame
* T: kinetic energy in a given reference frame
* equations: returns a vector of (symbolic) equations, m*a = F, projected
  along the vectors of one basis
* inertia: tensor of inertia about a point
* subs: takes values of the symbolic unknowns and returns a particle
  object which is purely numeric

AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef particle < anakin.point
    properties (Hidden = true, Access = protected)  
        M anakin.tensor = anakin.tensor(1); % mass of the object
    end
    properties
        forces cell = {}; % cell array with all forces (vectors) acting on the object
    end
    methods % creation
        function P = particle(varargin) % constructor
            switch nargin
                case 0 % no arguments
                    return;
                case 1 % particle                  
                    P.M = varargin{1}.M; 
                    P.v = varargin{1}.v; 
                case 2 % mass, vector 
                    P.M = anakin.tensor(varargin{1}); 
                    P.v = anakin.tensor(varargin{2});
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in particle');
            end
            P.curr_vel0 = P.vel;       
        end  
        function P = set.M(P,value) % on setting M
            P.M = anakin.tensor(value); 
        end
        function P = set.forces(P,value) % on setting forces
            if ~iscell(value) % ensure cell array
                value = {value};
            end
            for j=1:length(value) % validate input
                if ~isa(value{j},'anakin.tensor') 
                    error('The forces must be supplied in a cell array of anakin.tensor vectors');
                end
            end
            P.forces = value; 
        end
    end 
    methods (Hidden = true) % overloads
        function value = eq(P1,P2) % overload ==
            value = (P1.M == P2.M) && (P1.v == P2.v);
        end
        function value = ne(P1,P2) % overload =~
            value = ~eq(P1,P2);
        end
        function disp(P) % display
            disp('Point particle with mass:')
            disp(P.M.components)            
            disp('And coordinates:')
            disp(P.coordinates)            
        end
    end
    methods % general functionality 
        function mass = mass(P) % mass of the particle
            mass = P.M;
        end
        function p = p(P,S1) % linear momentum in S1
            if exist('S1','var')
                p = P.M*P.vel(S1);
            else
                p = P.M*P.vel;
            end            
        end
        function H = H(P,O,S1) % angular momentum about a point O in S1
            if ~exist('O','var')
                O = anakin.point; % default point
            end
            if exist('S1','var')
                H = cross(P.pos-O.pos, P.p(S1));
            else
                H = cross(P.pos-O.pos, P.p);
            end            
        end
        function T = T(P,S1) % kinetic energy in S1
            if exist('S1','var')
                vel = P.vel(S1).components;
            else
                vel = P.vel.components;
            end
            T = (P.M/2) * norm(vel)^2; 
        end
        function inertia = inertia(P,O) % inertia tensor of the particle with respect to point O
            if exist('O','var')
                r = P.pos - O.pos;
                inertia = P.M * (norm(r)^2 * anakin.tensor(eye(3)) - product(r,r));
            else
                inertia = anakin.tensor([0,0,0;0,0,0;0,0,0]); % no inertia about the particle itself
            end            
        end
        function eqs = equations(P,B1) % returns vector of equations of motion projected in basis B1
            p = P.p.dt;
            F = anakin.tensor([0;0;0]); % allocate;
            for i=1:length(P.forces)
                F = F + P.forces{i};
            end
            if ~exist('B1','var')
                B1 = anakin.basis;
            elseif isa(B1,'anakin.frame') 
                B1 = B1.basis; % extract basis
            end
            eqs = sym([0;0;0]); % allocate
            for i=1:3
                p_ = p * B1.e(i);
                F_ = F * B1.e(i);
                eqs(i) = (p_.components == F_.components);
            end                
        end
        function P_ = subs(P,variables,values) % particularize symbolic particle
            P_ = P;
            P_.v = P.v.subs(variables,values);
            P_.M = P.M.subs(variables,values);
            P_.curr_vel0 = P.curr_vel0.subs(variables,values);
        end
    end      
end




