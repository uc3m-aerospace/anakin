%{
DESCRIPTION:
particle: class to model a point particle. 

SYNTAX:
P0 = anakin.particle();  % returns default object  
P  = anakin.particle(<P>|(<mass>,<A|a|c>),<S1>);
where:
- <> denotes optional arguments
- | denotes alternative arguments
- P0 is the default particle (mass 1, located at the origin)
- P  is a particle  
- mass is a scalar or number
- A is a point
- a is a vector 
- c is an array with the Cartesian coordinates of the origin
- S1 is a frame. If given, all previous input as relative to that frame
 
METHODS: 
* mass: the mass of the particle
* p: linear momentum in a given reference frame
* H: angular momentum about a point in a given reference frame
* T: kinetic energy in a given reference frame
* I: tensor of inertia about a point
* subs: takes values of the symbolic unknowns and returns a particle
  object which is purely numeric

AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef particle < anakin.point
    properties (Hidden = true, Access = protected)  
        M anakin.tensor = anakin.tensor(1); % mass of the object
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
            b.v = S1.v; 
            for i = 1:length(varargin) % later inputs overwrite former inputs
                temp = varargin{i};
                if isa(temp,'anakin.body')
                    b.M = temp.M;                    
                    b.v = anakin.tensor(S1.v.components + temp.v.components(S1));  
                elseif isa(temp,'anakin.particle')
                    b.M = temp.M;
                    b.v = anakin.tensor(S1.v.components + temp.v.components(S1));                    
                elseif isa(temp,'anakin.point') % includes frame as a subclass
                    b.v = anakin.tensor(S1.v.components + temp.v.components(S1)); 
                elseif isa(temp,'anakin.tensor')
                    if temp.ndims == 0 % assume it is M
                        b.M = temp;
                    elseif temp.ndims == 1 % assume it is v
                        b.v = anakin.tensor(S1.v.components + temp.components(S1));
                    else
                        error('Cannot take tensors of order higher than 1 as inputs');
                    end
                else % Array
                    if numel(temp) == 1 % assume it is M
                        b.M = temp;
                    else % assume it is v
                        v_ = anakin.tensor(temp);
                        b.v = anakin.tensor(S1.v.components + v_.components(S1));
                    end
                end
            end 
        end
        function P = set.M(P,value) % on setting M
            P.M = anakin.tensor(value); 
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
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end 
            H = cross(P.pos-O.pos, P.p(S1));         
        end
        function T = T(P,S1) % kinetic energy in S1
            if ~exist('S1','var')
                S1 = anakin.frame; % default frame
            end 
            T = (P.M/2) * norm(P.vel(S1))^2; 
        end
        function I = I(P,O) % inertia tensor of the particle with respect to point O in canonical vector basis
            if ~exist('O','var')
                O = anakin.point; % default point is the origin
            end
            r = P.v - O.v;
            I = P.M * (norm(r)^2 * eye(3) - product(r,r));              
        end
        function P_ = subs(P,variables,values) % particularize symbolic particle
            P_ = P;
            P_.v = P.v.subs(variables,values);
            P_.M = P.M.subs(variables,values);
        end
    end      
end




