%{
DESCRIPTION:
particle: class to model a point particle. Inherits from point.

SYNTAX:
P0 = anakin.particle();  % returns default object  
P  = anakin.particle(<mass>,<P|A|a|c|>,<S1>);
where:
- <> denotes optional arguments
- | denotes alternative arguments
- P0 is the default particle (mass 1, located at the origin)
- P  is a particle 
- mass is the object mass
- A is a point
- a is a vector
- c is an array with the three vector components 
- S1 is a frame. If given, all previous input as relative to that frame

PROPERTIES:
* mass

METHODS: 
* p: linear momentum in a given reference frame
* H: angular momentum about a point in a given reference frame
* T: kinetic energy in a given reference frame
* inertia: tensor of inertia about a point
* subs: takes values of the symbolic unknowns and returns a particle
  object which is purely numeric

AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef particle < anakin.point 
    properties  
        mass anakin.tensor = anakin.tensor(1); % mass of the object
    end
    methods % creation
        function P = particle(varargin) % constructor 
            switch nargin
                case 0 % no arguments
                    return;
                case 1  
                    Pt = anakin.particle(varargin{1},anakin.frame);
                    P.pos0 = Pt.pos0; 
                    P.mass = Pt.mass;
                case 2  
                    if isa(varargin{end},'anakin.frame') % last is frame
                        if isa(varargin{1},'anakin.point') || isa(varargin{1},'anakin.tensor') || numel(varargin{1}) == 3 % (vector or components), frame
                            P.pos0 = anakin.tensor(varargin{2}.basis.matrix * anakin.tensor(varargin{1}).components + varargin{2}.origin.coordinates);
                        else % mass, frame
                            P.pos0 = varargin{2}.origin;
                            P.mass = anakin.tensor(varargin{1});
                        end
                    else
                        Pt = anakin.particle(varargin{1},varargin{2},anakin.frame);
                        P.pos0 = Pt.pos0;
                        P.mass = Pt.mass;
                    end 
                case 3   
                    P.pos0 = anakin.tensor(varargin{3}.basis.matrix * anakin.tensor(varargin{2}).components + varargin{3}.origin.coordinates);
                    P.mass = anakin.tensor(varargin{1});
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in particle');
            end       
        end 
        function P = set.mass(P,value) % on setting c
            P.mass = value;
            if isa(P.mass,'sym') % symbolic input
                P.mass = formula(simplify(P.mass)); % simplify and force sym rather than symfun to allow indexing
            end
            try
                P.mass = double(P.mass);
            end
        end
    end 
    methods
        function disp(P) % display
            disp('Particle with mass:')
            disp(P.mass.components)            
            disp('canonical position vector components:')
            disp(P.pos.components)            
        end
    end
    methods % functionality 
        function p = p(P,S1) % linear momentum in S1
            if exist('S1','var')
                p = P.mass*P.vel(S1);
            else
                p = P.mass*P.vel;
            end            
        end
        function H = H(P,O,S1) % angular momentum about O in S1
            if ~exist('O','var')
                O = anakin.point; % default point
            end
            if exist('S1','var')
                H = cross(P.pos0-O.pos0, P.p(S1));
            else
                H = cross(P.pos0-O.pos0, P.p);
            end            
        end
        function T = T(P,S1) % kinetic energy in S1
            if exist('S1','var')
                vel = P.vel(S1).components;
            else
                vel = P.vel.components;
            end
            T = (P.mass/2) * dot(vel,vel); 
        end
        function inertia = inertia(P,O) % inertia tensor of the particle with respect to point O
            if exist('O','var')
                r = P.pos0 - O;
            else
                r = P.pos0;
            end
            inertia = P.mass * (norm(r)^2 * anakin.tensor(eye(3)) - product(r,r));
        end
        function P_ = subs(P,variables,values) % particularize symbolic vector
            P_ = P;
            P_.pos0 = P.pos0.subs(variables,values);
            P_.mass = P.mass.subs(variables,values);
        end
    end      
end




