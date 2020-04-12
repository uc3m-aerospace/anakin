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
        function S = frame(varargin) % constructor
            switch nargin
                case 0 % no arguments
                    return;
                case 1
                    S = anakin.frame(varargin{1},anakin.frame);
                case 2
                    if isa(varargin{end},'anakin.frame') 
                        if isa(varargin{1},'anakin.frame')
                            S.v = anakin.tensor(varargin{2}.v.components + varargin{2}.m * varargin{1}.v.components); 
                            S.m = varargin{2}.m * varargin{1}.m; 
                        elseif isa(varargin{1},'anakin.point')
                            S.v = anakin.tensor(varargin{2}.v.components + varargin{2}.m * varargin{1}.v.components); 
                        elseif isa(varargin{1},'anakin.basis')
                            S.m = varargin{2}.m * varargin{1}.m;        
                        elseif isa(varargin{1},'anakin.tensor')
                            S.v = anakin.tensor(varargin{2}.v.components + varargin{2}.m * varargin{1}.components); 
                        elseif length(varargin{1}(:,1)) > 1 && length(varargin{1}(1,:)) > 1
                            S.m = varargin{2}.m * varargin{1}; 
                        else
                            temp = anakin.tensor(varargin{1});
                            S.v = anakin.tensor(varargin{2}.v.components + varargin{2}.m * temp.components); 
                        end
                    else
                        S = anakin.frame(varargin{1},varargin{2},anakin.frame);
                    end
                case 3
                    if isa(varargin{1},'anakin.point')
                        S.v = anakin.tensor(varargin{3}.v.components + varargin{3}.m * varargin{1}.v.components);
                    elseif isa(varargin{1},'anakin.tensor')
                        S.v = anakin.tensor(varargin{3}.v.components + varargin{3}.m * varargin{1}.components);
                    else
                        temp = anakin.tensor(varargin{1});
                        S.v = anakin.tensor(varargin{3}.v.components + varargin{3}.m * temp.components);
                    end
                    if isa(varargin{2},'anakin.basis')
                        S.m = varargin{3}.m * varargin{2}.m; 
                    else 
                        S.m = varargin{3}.m * varargin{2};
                    end
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in frame');
            end     
        end      
    end
    methods (Hidden = true) % overloads
        function value = eq(S1,S2) % overload ==
            value = (S1.v == S2.v) && (S1.m == S2.m);
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
            S_.b = S.b.subs(variables,values);
        end         
    end 
    methods % 3-dimensional-space functionality 
        function h = plot(S,varargin) % plot   
            h = [S.origin.plot(varargin{:});S.basis.plot(S.origin,varargin{:})];
        end
    end 
end
