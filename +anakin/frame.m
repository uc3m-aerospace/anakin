%{
DESCRIPTION:
frame: class to define reference frames

SYNTAX:
S0 = anakin.frame();  % returns default object 
S  = anakin.frame(<S|(<A|a|c>,<B|m>)>,<S1>);
where:
- <> denotes optional arguments
- | denotes alternative arguments
- () groups argument options
- S0 is the default frame (canonical reference frame)
- S is a frame
- A is a point
- a is a vector (1st-order tensor) 
- c is an array with the Cartesian coordinates of the origin
- B  is a basis 
- m  is a square matrix (square array) for the basis
- S1 is a frame. If given, all previous input as relative to that frame
  
METHODS: 
* origin, basis: these methods return the origin and the basis of the frame
* subs: takes values of the symbolic unknowns and returns a reference frame with
  purely numeric origin point and basis matrix (symbolic variables must be used)     
* plot: plots the reference frame

AUTHOR:
Mario Merino <mario.merino@uc3m.es>
%}
classdef frame < anakin.point & anakin.basis 
    methods % creation 
        function b = frame(varargin) % constructor
            if isempty(varargin) % Default
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
                if isa(temp,'anakin.frame') % includes body as subclass
                    b.v = anakin.tensor(S1.v.components + temp.v.components(S1));
                    b.m = temp.matrix(S1);
                elseif isa(temp,'anakin.point') % includes particle as subclass
                    b.v = anakin.tensor(S1.v.components + temp.v.components(S1)); 
                elseif isa(temp,'anakin.basis')
                    b.m = temp.matrix(S1);
                elseif isa(temp,'anakin.tensor')
                    b.v = anakin.tensor(S1.v.components + temp.components(S1));
                else % Array 
                    v_ = anakin.tensor(temp);
                    b.v = anakin.tensor(S1.v.components + v_.components(S1)); 
                end
            end 
        end
    end
    methods (Hidden = true) % overloads
        function value = eq(S1,S2) % overload ==
            value = (S1.v == S2.v) && (S1.basis == S2.basis);
        end
        function value = ne(S1,S2) % overload ~=
            value = ~eq(S1,S2);
        end
        function disp(S) % display
            disp('Frame with origin with coordinates:')
            disp(S.coordinates)
            disp('And basis with rotation matrix:')
            disp(S.matrix)
        end
    end
    methods % general functionality    
        function origin = origin(S) % return the origin point
            origin = anakin.point(S.v);
        end
        function basis = basis(S) % return the basis
            basis = anakin.basis(S.m);
        end 
        function S_ = subs(S,variables,values) % Particularize symbolic frame
            S_ = S;   
            S_.v = S.v.subs(variables,values); 
            S_.m = subs(S.m,variables,values);
            try
                S_.m = double(S_.m);
            catch
                % pass
            end
            
        end         
    end 
    methods % 3-dimensional-space functionality 
        function h = plot(S,varargin) % plot   
            h = [S.origin.plot(varargin{:});S.basis.plot(S.origin,varargin{:})];
        end
    end 
end
