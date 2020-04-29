%{
DESCRIPTION:
frame: class to define reference frames. Subclass of basis and point.

SYNTAX:
S = anakin.frame();  % returns default object (canonical reference frame)
S = anakin.frame(...,<S1>); 
where:
- <> denotes optional arguments 
- S is a frame
- ... denotes a list of one or more of the following. Later inputs can
  overwrite previous inputs:
    - frame object
    - point object: denotes the origin point of the reference frame
    - vector (1st-order tensor): denotes the origin point of the reference
      frame 
    - one dimensional array: denotes the origin point of the reference
      frame 
    - basis: denotes the basis of the reference frame
    - A square rotation matrix: denotes the basis of the reference frame
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
            b.r = S1.r;
            b.m = S1.m;
            for i = 1:length(varargin) % later inputs overwrite former inputs
                temp = varargin{i}; 
                if isa(temp,'anakin.frame') % includes body as subclass
                    b.r = anakin.tensor(S1.r.components + temp.r.components(S1));
                    b.m = temp.matrix(S1);
                elseif isa(temp,'anakin.point') % includes particle as subclass
                    b.r = anakin.tensor(S1.r.components + temp.r.components(S1)); 
                elseif isa(temp,'anakin.basis')
                    b.m = temp.matrix(S1);
                elseif isa(temp,'anakin.tensor')
                    b.r = anakin.tensor(S1.r.components + temp.components(S1));
                else % Array  
                    v_ = anakin.tensor(temp);
                    if v_.ndims == 1 % it is a vector
                        b.r = anakin.tensor(S1.r.components + v_.components(S1)); 
                    else % it is a basis
                        b.m = anakin.basis(v_.components,S1).matrix;
                    end
                end
            end 
        end
    end
    methods (Hidden = true) % overloads
        function value = eq(S1,S2) % overload ==
            value = (S1.r == S2.r) && (S1.basis == S2.basis);
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
            origin = anakin.point(S.r);
        end
        function basis = basis(S) % return the basis
            basis = anakin.basis(S.m);
        end 
        function S_ = subs(S,variables,values) % Particularize symbolic frame
            S_ = S;   
            S_.r = S.r.subs(variables,values); 
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
