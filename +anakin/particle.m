%{
DESCRIPTION:
particle: class to model a point particle. Subclass of point.

SYNTAX:
P = anakin.particle();  % returns default object (mass 1, at the origin) 
P = anakin.particle(...,<S1>);
where:
- <> denotes optional arguments
- P is a particle
- ... denotes a list of one or more of the following. Later inputs can
  overwrite previous inputs: 
    - particle object
    - point object: denotes the origin point of the reference frame
    - scalar (0th-order tensor): mass
    - scalar numeric or symbolic value: mass
    - vector (1st-order tensor): denotes the position  
    - one dimensional array: denotes the position  
- S1 is a frame. If given, all previous inputs are relative to that frame
 
PROPERTIES:
* mass: the mass of the particle

METHODS: 

* p: linear momentum in a given reference frame
* H: angular momentum about a point in a given reference frame
* T: kinetic energy in a given reference frame
* I: tensor of inertia about a point
* subs: takes values of the symbolic unknowns and returns a particle
  object which is purely numeric
* force_equation: returns the force equation projected along a desired
  direction. Requires symbolic particle

AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef particle < anakin.point
    properties 
        mass anakin.tensor = anakin.tensor(1); % mass of the object 
    end 
    methods % creation 
        function b = particle(varargin) % constructor
            if isempty(varargin) % Default
                return;
            elseif length(varargin) > 1 && isa(varargin{end},'anakin.frame') % Last argin is frame
                S1 = varargin{end};
                varargin = varargin(1:end-1); 
            else % No frame is provided; use default
                S1 = anakin.frame;
            end 
            b.r = S1.r; 
            for i = 1:length(varargin) % later inputs overwrite former inputs
                temp = varargin{i};
                if isa(temp,'anakin.body')
                    b.mass = temp.mass;                    
                    b.r = anakin.tensor(S1.r.components + temp.r.components(S1));  
                elseif isa(temp,'anakin.particle')
                    b.mass = temp.mass;
                    b.r = anakin.tensor(S1.r.components + temp.r.components(S1));                    
                elseif isa(temp,'anakin.point') % includes frame as a subclass
                    b.r = anakin.tensor(S1.r.components + temp.r.components(S1)); 
                elseif isa(temp,'anakin.tensor')
                    if temp.ndims == 0 % assume it is mass
                        b.mass = temp;
                    elseif temp.ndims == 1 % assume it is r
                        b.r = anakin.tensor(S1.r.components + temp.components(S1));
                    else
                        error('Cannot take tensors of order higher than 1 as inputs');
                    end
                else % Array
                    v_ = anakin.tensor(temp);
                    if v_.ndims == 0 % assume it is mass
                        b.mass = v_;
                    elseif v_.ndims == 1 % assume it is r
                        b.r = anakin.tensor(S1.r.components + v_.components(S1));
                    else
                        error('Cannot take arrays of order higher than 1 as inputs');
                    end
                end
            end 
        end
        function P = set.mass(P,value) % on setting mass
            P.mass = anakin.tensor(value); 
            if P.mass.ndims ~= 0
                error('Mass must be a scalar');
            end
        end 
    end 
    methods (Hidden = true) % overloads
        function value = eq(P1,P2) % overload ==
            value = (P1.mass == P2.mass) && (P1.r == P2.r);
        end
        function value = ne(P1,P2) % overload =~
            value = ~eq(P1,P2);
        end
        function disp(P) % display
            disp('Point particle with mass:')
            disp(P.mass.components)            
            disp('And coordinates:')
            disp(P.coordinates)            
        end
    end
    methods % general functionality  
        function p = p(P,S1) % linear momentum in S1
            if exist('S1','var')
                p = P.mass*P.vel(S1);
            else
                p = P.mass*P.vel;
            end            
        end
        function H = H(P,O,S1) % angular momentum about a point O in S1
            if ~exist('O','var')
                O = anakin.point; % default point
            end
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end 
            H = cross(P.pos-O.pos, P.p(S1));         
        end
        function T = T(P,S1) % kinetic energy in S1
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end 
            T = (P.mass/2) * norm(P.vel(S1))^2; 
        end
        function I = I(P,O) % inertia tensor of the particle with respect to point O in canonical vector basis
            if ~exist('O','var')
                O = anakin.point; % default point is the origin
            end
            r = P.r - O.r;
            I = P.mass * (norm(r)^2 * eye(3) - product(r,r));              
        end
        function P_ = subs(P,variables,values) % particularize symbolic particle
            P_ = P;
            P_.r = P.r.subs(variables,values);
            P_.mass = P.mass.subs(variables,values);
        end
    end    
    methods % dynamics
        function eq = force_equation(P,F,e,S) % force equation along direction e, assuming S is inertial. Requires symbolic particle
            if ~exist('S','var')
                S = anakin.frame; % default to canonical frame
            end
            lhs = P.p(S).dt(S)*e;
            rhs = F*e; 
            eq = sym(lhs.components - rhs.components); % expression equal to zero is the equation
        end  
    end
end




