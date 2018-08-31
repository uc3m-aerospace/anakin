%{
particle: class to model a point particle. Inherits from point.

P0 = anakin.particle();  % returns default object  
P  = anakin.particle(<mass>,<P|A|a|c|x,y,z>,<S1>);

where:
- <> denotes optional arguments
- | denotes alternative arguments
- P0 is the default particle (mass 1, located at the origin)
- P  is a particle 
- mass is the object mass
- A is a point
- a is a vector
- c is an array with the three vector components
- x,y,z are the three vector components
- S1 is a frame. If given, all previous input as relative to that frame

METHODS: 
* p: linear momentum in a given reference frame
* H: angular momentum about a point in a given reference frame
* T: kinetic energy in a given reference frame
* inertia: tensor of inertia about a point
* subs: takes values of the symbolic unknowns and returns a particle
  object which is purely numeric

AUTHOR: Mario Merino <mario.merino@uc3m.es>
%}
classdef particle < anakin.point 
    properties (SetAccess = protected)
        mass = 1; % mass of the object
    end
    methods % creation
        function P = particle(varargin) % constructor
            for i = 1:length(varargin)
               if isa(varargin{i},'sym')
                   varargin{i} = formula(varargin{i}); % enforce formula to allow indexing
               end
            end            
            switch nargin
                case 0 % no arguments
                    return;
                case 1  
                    Pt = anakin.particle(varargin{1},anakin.frame);
                    P.c = Pt.c; 
                    P.mass = Pt.mass;
                case 2  
                    if isa(varargin{end},'anakin.frame') % last is frame
                        if isa(varargin{1},'anakin.vector') || numel(varargin{1}) == 3 % (vector or components), frame
                            P.c = varargin{2}.matrix * anakin.vector(varargin{1}).c + varargin{2}.c;
                        else % mass, frame
                            P.c = varargin{2}.c;
                            P.mass = varargin{1};
                        end
                    else
                        Pt = anakin.particle(varargin{1},varargin{2},anakin.frame);
                        P.c = Pt.c;
                        P.mass = Pt.mass;
                    end 
                case 3  
                    if isa(varargin{end},'anakin.frame') % last is frame
                        if isa(varargin{2},'anakin.vector') || numel(varargin{2}) == 3 % (vector or components), mass, frame
                            P.c = varargin{3}.matrix * anakin.vector(varargin{2}).c + varargin{3}.c;
                            P.mass = varargin{1}; 
                        end
                    else
                        Pt = anakin.particle(varargin{1},varargin{2},varargin{3},anakin.frame);
                        P.c = Pt.c;
                        P.mass = Pt.mass;
                    end  
                case 4 
                    if isa(varargin{end},'anakin.frame') % last is frame 
                        P.c = varargin{4}.matrix * anakin.vector(varargin{1},varargin{2},varargin{3}).c + varargin{4}.c; 
                    else
                        Pt = anakin.particle(varargin{1},varargin{2},varargin{3},varargin{4},anakin.frame);
                        P.c = Pt.c;
                        P.mass = Pt.mass;
                    end   
                case 5 % xyz, mass, frame
                    P.c = varargin{5}.matrix * anakin.vector(varargin{2},varargin{3},varargin{4}).c + varargin{5}.c;
                    P.mass = varargin{1}; 
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in particle');
            end       
        end 
        function P = set.mass(P,value) % on setting c
            P.mass = value;
            if isa(P.mass,'sym') % symbolic input
                P.mass = formula(simplify(P.mass)); % simplify and force sym rather than symfun to allow indexing
            end
        end
    end 
    methods
        function disp(P) % display
            disp('Point particle with mass:')
            disp(P.mass)            
            disp('and canonical position:')
            disp(P.c)            
        end
    end
    methods % functionality 
        function p = p(P,S1) % linear momentum in S1
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end
            p = P.mass*P.vel(S1);
        end
        function H = H(P,O,S1) % angular momentum about O in S1
            if ~exist('O','var')
                O = anakin.point; % default point
            end
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end
            H = cross(anakin.vector(P)-anakin.vector(O), P.mass*P.vel(S1));
        end
        function T = T(P,S1) % kinetic energy in S1
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end
            T = (P.mass/2) * dot(P.vel(S1),P.vel(S1));
            if isa(T,'sym') % symbolic input
                T = formula(simplify(T)); % simplify and force sym rather than symfun to allow indexing
            end
        end
        function P_ = subs(P,variables,values) % particularize symbolic vector
            P_ = P;
            P_.c = double(subs(P.c,variables,values));
            P_.mass = double(subs(P.mass,variables,values));
        end
    end      
end




