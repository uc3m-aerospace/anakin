%{
DESCRIPTION:
frame: class to define reference frames

SYNTAX:
S0 = anakin.frame();  % returns default object
S  = anakin.frame(S,<S1>); % (convert to class)
S  = anakin.frame(<A|a|c>,<B|m>,<S1>);
where:
- <> denotes optional arguments
- | denotes alternative arguments
- () groups argument options
- S0 is the default frame (canonical reference frame)
- A is a point
- a is a vector
- c is an array with the three vector components
- B  is a basis
- m  is a matrix
- S1 is a frame. If given, all previous input as relative to that frame
 
PROPERTIES:
* origin: origin point (anakin.tensor)
* basis: reference frame basis (anakin.basis)

METHODS:
* spacedim: returns dimensionality of space
* subs: takes values of the symbolic unknowns and returns a reference frame with
  purely numeric origin point and basis matrix (symbolic variables must be used)     
* plot: plots the frame with quiver

AUTHOR:
Mario Merino <mario.merino@uc3m.es>
%}
classdef frame
    properties
        origin anakin.point = anakin.point;
        basis  anakin.basis = anakin.basis;        
    end
    methods % creation
        function S = frame(varargin) % constructor 
            switch nargin
                case 0 % no arguments
                    return;
                case 1
                    St = anakin.frame(varargin{1},anakin.frame);
                    S.origin = St.origin;
                    S.basis = St.basis;                     
                case 2
                    if isa(varargin{end},'anakin.frame') % last is frame
                        if isa(varargin{1},'anakin.frame') % relative frame, frame
                            S.origin = anakin.tensor(varargin{2}.basis.matrix * varargin{1}.origin.coordinates + varargin{2}.origin.coordinates); 
                            S.basis = anakin.basis(varargin{2}.basis.matrix * varargin{1}.basis.matrix);
                        elseif isa(varargin{1},'anakin.point') || isa(varargin{1},'anakin.tensor') || numel(varargin{1}) == 3 % (relative vector or relative column), frame
                            S.origin = anakin.tensor(varargin{2}.basis.matrix * anakin.tensor(varargin{1}).components + varargin{2}.origin.coordinates); 
                            S.basis = varargin{2}.basis; % copy basis from given frame 
                        else % (relative basis or relative matrix), frame
                            S.origin = varargin{2}.origin; % copy origin from given frame
                            S.basis = anakin.basis(varargin{2}.basis.matrix * anakin.basis(varargin{1}).matrix);  
                        end
                    else 
                        St = anakin.frame(varargin{1},varargin{2},anakin.frame);
                        S.origin = St.origin;
                        S.basis = St.basis;      
                    end
                case 3 % (relative vector or relative column), (relative basis or relative matrix), frame                         
                    S.origin = anakin.tensor(varargin{3}.basis.matrix * anakin.tensor(varargin{1}).components + varargin{3}.origin.coordinates);
                    S.basis = anakin.basis(varargin{3}.basis.matrix * anakin.basis(varargin{2}).matrix);
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in frame');
            end     
        end        
    end
    methods (Hidden = true) % overloads
        function value = eq(S1,S2) % overload ==
            value = (S1.origin == S2.origin) && (S1.basis == S2.basis);
        end
        function value = ne(S1,S2) % overload ~=
            value = ~eq(S1,S2);
        end
        function disp(S) % display
            disp('Frame with origin with canonical position:')
            disp(S.origin.components)
            disp('and basis with canonical rotation matrix:')
            disp(S.basis.matrix)
        end
    end
    methods % general functionality   
        function spacedim = spacedim(S) % number of dimensions of space
            spacedim = S.basis.spacedim;
            if spacedim ~= S.origin.spacedim
                error('Inconsistent space dimensionality in the origin and basis of the reference frame')
            end
        end            
        function S_ = subs(S,variables,values) % Particularize symbolic frame
            S_ = S;   
            S_.origin = S.origin.subs(variables,values);
            S_.basis = S.basis.subs(variables,values);
        end         
    end 
    methods % 3-dimensional-space functionality 
        function h = plot(S,varargin) % plot  
            h = [S.origin.plot(varargin{:}), S.basis.plot(S.origin.pos,varargin{:})];
        end
    end 
end
