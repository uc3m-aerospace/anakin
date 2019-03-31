%{
DESCRIPTION:
frame: class to define reference frames

SYNTAX:
S0 = anakin.frame();  % returns default object 
S  = anakin.frame(<S|(<A|a>,<B>)>,<S1>);
where:
- <> denotes optional arguments
- | denotes alternative arguments
- () groups argument options
- S0 is the default frame (canonical reference frame)
- S is a frame
- A is a point
- a is a vector (1st-order tensor) 
- B  is a basis 
- S1 is a frame. If given, all previous input as relative to that frame
 
PROPERTIES:
* origin: origin point (anakin.tensor)
* basis: reference frame basis (anakin.basis)

METHODS: 
* origin, basis: these methods return the origin and the basis of the frame
* subs: takes values of the symbolic unknowns and returns a reference frame with
  purely numeric origin point and basis matrix (symbolic variables must be used)     
* plot: plots the frame with quiver

AUTHOR:
Mario Merino <mario.merino@uc3m.es>
%}
classdef frame < anakin.point
    properties (Hidden = true, Access = protected) 
        b anakin.basis = anakin.basis;        
    end
    methods % creation
        function S = frame(varargin) % constructor
            switch nargin
                case 0 % no arguments
                    return;
                case 1 % frame
                    S.v = varargin{1}.v;
                    S.b = varargin{1}.b;
                case 2 % vector,
                    S.v = anakin.tensor(varargin{1});
                    S.b = anakin.basis(varargin{2});                     
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in frame');
            end     
        end
        function S = set.b(S,value) % on setting b
            S.b = anakin.basis(value); 
        end        
    end
    methods (Hidden = true) % overloads
        function value = eq(S1,S2) % overload ==
            value = (S1.v == S2.v) && (S1.b == S2.b);
        end
        function value = ne(S1,S2) % overload ~=
            value = ~eq(S1,S2);
        end
        function disp(S) % display
            disp('Frame with origin with coordinates:')
            disp(S.coordinates)
            disp('And basis with rotation matrix:')
            disp(S.b.matrix)
        end
    end
    methods % general functionality    
        function origin = origin(S) % return the origin point
            origin = anakin.point(S.v);
        end
        function basis = basis(S) % return the basis
            basis = S.b; 
        end
        function omega = omega(S,S1) % omega vector with respect to reference frame S1
            omega = S.b.omega(S1.b);
        end
        function alpha = alpha(S,S1) % alpha vector with respect to reference frame S1
            alpha = S.b.alpha(S1.b);
        end 
        function S_ = subs(S,variables,values) % Particularize symbolic frame
            S_ = S;   
            S_.v = S.v.subs(variables,values);
            S_.b = S.b.subs(variables,values);
        end         
    end 
    methods % 3-dimensional-space functionality 
        function h = plot(S,varargin) % plot  
            h = [S.origin.plot(varargin{:}), S.b.plot(S.origin.pos,varargin{:})];
        end
    end 
end
