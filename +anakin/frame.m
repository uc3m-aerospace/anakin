%{
frame: class to define orthogonal, right-handed reference frames.
Inherits from point and basis.

S0 = anakin.frame();  % no arguments return default object
S  = anakin.frame(S,<S1>); % (convert to class)
S  = anakin.frame(<A|a|c|x,y,z>,<B|m|(a|c|x,y,z),(a|c|x,y,z),(a|c|x,y,z)|q|axis,angle>,<S1>);

where:
- <> denotes optional arguments
- | denotes alternative arguments
- () groups argument options
- S0 is the default frame (canonical reference frame)
- A is a point
- a is a vector
- c is an array with the three vector components
- x,y,z are the three vector components
- B  is a basis
- m  is a matrix
- q are quaternions
- axis is the unit vector of the axis of rotation
- angle is angle of rotation about axis
- S1 is a frame. If given, all previous input as relative to that frame
 
METHODS:
* origin: returns the origin point
* basis: returns the basis
* subs: takes values of the symbolic unknowns and returns a reference frame with
  purely numeric origin point and basis matrix (symbolic variables must be used)     
* plot: plots the frame with quiver

AUTHOR: Mario Merino <mario.merino@uc3m.es>
%}
classdef frame < anakin.point & anakin.basis % Inherit from point and basis 
    methods % creation
        function S = frame(varargin) % constructor
            for i = 1:length(varargin)
               if isa(varargin{i},'sym')
                   varargin{i} = formula(varargin{i}); % enforce formula to allow indexing
               end
            end
            switch nargin
                case 0 % no arguments
                    return;
                case 1
                    St = anakin.frame(varargin{1},anakin.frame);
                    S.c = St.c;
                    S.m = St.m;                     
                case 2
                    if isa(varargin{end},'anakin.frame') % last is frame
                        if isa(varargin{1},'anakin.frame') % relative frame, frame
                            S.c = varargin{2}.m * varargin{1}.c + varargin{2}.c; 
                            S.m = varargin{2}.m * varargin{1}.m; % copy basis from given frame 
                        elseif isa(varargin{1},'anakin.vector') || numel(varargin{1}) == 3 % (relative vector or relative column), frame
                            S.c = varargin{2}.m * anakin.vector(varargin{1}).c + varargin{2}.c; 
                            S.m = varargin{2}.m; % copy basis from given frame 
                        else % (relative basis or relative matrix), frame
                            S.c = varargin{2}.c; % copy origin from given frame
                            S.m = varargin{2}.m * anakin.basis(varargin{1}).m;  
                        end
                    else 
                        St = anakin.frame(varargin{1},varargin{2},anakin.frame);
                        S.c = St.c;
                        S.m = St.m;
                    end
                case 3 
                    if isa(varargin{end},'anakin.frame') % (relative vector or relative column), (relative basis or relative matrix), frame 
                        if isa(varargin{2},'anakin.basis') || numel(varargin{2}) == 9 || numel(varargin{2}) == 4 % (vector or column), (basis, matrix or quaternions), frame
                            S.c = varargin{3}.m * anakin.vector(varargin{1}).c + varargin{3}.c;
                            S.m = varargin{3}.m * anakin.basis(varargin{2}).m;  
                        else % axis, angle, frame
                            S.c = varargin{3}.c;
                            S.m = varargin{3}.m * anakin.basis(varargin{1},varargin{2}).m;  
                        end
                    else
                        St = anakin.frame(varargin{1},varargin{2},varargin{3},anakin.frame);
                        S.c = St.c;
                        S.m = St.m;
                    end
                case 4  
                    if isa(varargin{end},'anakin.frame') % last is frame
                        if isa(varargin{3},'anakin.vector') || numel(varargin{3}) == 3 % relative ijk (vectors or columns), frame                        
                            S.c = varargin{4}.c; % copy origin from given frame
                            S.m = varargin{4}.m * anakin.basis(varargin{1},varargin{2},varargin{3}).m;                    
                        elseif isa(varargin{2},'anakin.vector') || numel(varargin{2}) == 3 % (vector or column), axis, angle, frame                        
                            S.c = varargin{4}.m * anakin.vector(varargin{1}).c + varargin{4}.c; 
                            S.m = varargin{4}.m * anakin.basis(varargin{2},varargin{3}).m;                                            
                        else % relative xyz, frame
                            S.c = varargin{4}.m * anakin.vector(varargin{1},varargin{2},varargin{3}).c + varargin{4}.c; 
                            S.m = varargin{4}.m; % copy basis from given frame
                        end
                    else
                        St = anakin.frame(varargin{1},varargin{2},varargin{3},varargin{4},anakin.frame);
                        S.c = St.c;
                        S.m = St.m;
                    end
                case 5
                    if isa(varargin{end},'anakin.frame') % last is frame
                        if isa(varargin{4},'anakin.basis') || numel(varargin{4}) == 9 || numel(varargin{4}) == 4 % relative xyz, (relative basis or relative matrix or quaternions), frame
                            S.c = varargin{5}.m * anakin.vector(varargin{1},varargin{2},varargin{3}).c + varargin{5}.c;
                            S.m = varargin{5}.m * anakin.basis(varargin{4}).m;
                        else % (relative vector or relative column), relative ijk, frame
                            S.c = varargin{5}.m * anakin.vector(varargin{1}).c + varargin{5}.c;
                            S.m = varargin{5}.m * anakin.basis(varargin{2},varargin{3},varargin{4}).m; 
                        end
                    else
                        St = anakin.frame(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},anakin.frame);
                        S.c = St.c;
                        S.m = St.m;
                    end
                case 6
                    if isa(varargin{end},'anakin.frame') % last is frame                        
                        S.c = varargin{6}.m * anakin.vector(varargin{1},varargin{2},varargin{3}).c + varargin{6}.c; % xyz, axis, angle, frame
                        S.m = varargin{6}.m * anakin.basis(varargin{4},varargin{5}).m;                         
                    else
                        St = anakin.frame(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},anakin.frame);
                        S.c = St.c;
                        S.m = St.m;
                    end
                case 7 % relative x,y,z, relative i,j,k, frame
                    S.c = varargin{7}.m * anakin.vector(varargin{1},varargin{2},varargin{3}).c + varargin{7}.c;
                    S.m = varargin{7}.m * anakin.basis(varargin{4},varargin{5},varargin{6}).m;
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in frame');
            end     
        end        
    end
    methods % overloads
        function value = eq(S1,S2) % overload ==
            value = (eq@anakin.point(S1,S2) && eq@anakin.basis(S1,S2));
        end
        function value = ne(S1,S2) % overload ~=
            value = ~eq(S1,S2);
        end
    end
    methods % functionality 
        function origin = origin(S) % returns the origin point
            origin = anakin.point(S);
        end
        function basis = basis(S) % returns the origin point
            basis = anakin.basis(S);
        end
        function S_ = subs(S,variables,values) % Particularize symbolic frame
            S_ = S;   
            S_.c = double(subs(S.c,variables,values));
            S_.m = double(subs(S.m,variables,values));
        end         
    end
    methods % plotting
        function h = plot(S,varargin) % plot  
            cc = S.c;
            mm = S.m;
            hold on            
            h = quiver3([cc(1),cc(1),cc(1)],[cc(2),cc(2),cc(2)],[cc(3),cc(3),cc(3)],...
            mm(1,:),mm(2,:),mm(3,:),0,'color','k');
            hold off
            if ~isempty(varargin)
                set(h,varargin{:}); % set options stored in varargin
            end
        end
    end
    methods % removed methods 
        function mtimes(~,~)
            error('This usage of point is not permitted');
        end 
        function mrdivide(~,~)
            error('This usage of point is not permitted');
        end 
        function mldivide(~,~)
            error('This usage of point is not permitted');
        end 
    end
end
