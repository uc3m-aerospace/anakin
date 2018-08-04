%{
frame: class to define reference frames

A frame can be created by passing a origin point and a vector basis with respect
to the canonical reference frame S0. Alternatively, a relative origin point
poisition, a relative basis vector basis and a different reference frame can be
given. Called without arguments, frame returns the canonical reference frame S0.   
This class inherits from point and basis.

Equality and non-equality operators have been overloaded to apply to
frames too.
 
METHODS:
* origin: returns the origin point
* subs: takes values of the symbolic unknowns and returns a reference frame with
  purely numeric origin point and basis matrix (symbolic variables must be used)     
* plot: plots the frame with quiver

MMM20180802
%}
classdef frame < anakin.point & anakin.basis % Inherit from point and basis 
    methods % creation
        function S = frame(varargin) % creator
            switch nargin
                case 0 % no arguments
                    return;
                case 1
                    if isa(varargin{1},'anakin.frame') % frame
                        S = varargin{1};
                    elseif isa(varargin{1},'anakin.vector') || numel(varargin{1}) == 3 % vector or column
                        S.c = anakin.vector(varargin{1}).c; 
                    else % basis or matrix
                        S.m = anakin.basis(varargin{1}).m;      
                    end
                case 2
                    if isa(varargin{2},'anakin.frame') % last is frame
                        if isa(varargin{1},'anakin.vector') || numel(varargin{1}) == 3 % (relative vector or relative column), frame
                            S.c = varargin{2}.m * anakin.vector(varargin{1}).c + varargin{2}.c; 
                            S.m = varargin{2}.m; % copy basis from given frame 
                        else % (relative basis or relative matrix), frame
                            S.c = varargin{2}.c; % copy origin from given frame
                            S.m = varargin{2}.m * anakin.basis(varargin{1}).m;  
                        end
                    else % (vector or column), (basis or matrix)
                        S.c = anakin.vector(varargin{1}).c;
                        S.m = anakin.basis(varargin{2}).m;
                    end
                case 3 
                    if isa(varargin{3},'anakin.frame') % (relative vector or relative column), (relative basis or relative matrix), frame 
                        S.c = varargin{3}.m * anakin.vector(varargin{1}).c + varargin{3}.c;
                        S.m = varargin{3}.m * anakin.basis(varargin{2}).m;                         
                    elseif isa(varargin{3},'anakin.vector') || numel(varargin{3}) == 3 % i,j,k
                        S.m = anakin.basis(varargin{1},varargin{2},varargin{3}).m;     
                    else % x,y,z
                        S.c = anakin.vector(varargin{1},varargin{2},varargin{3}).c; 
                    end
                case 4  
                    if isa(varargin{4},'anakin.frame') % last is frame
                        if isa(varargin{1},'anakin.vector') || numel(varargin{1}) == 3 % relative ijk (vectors or columns), frame                        
                            S.c = varargin{4}.c; % copy origin from given frame
                            S.m = varargin{4}.m * anakin.basis(varargin{1},varargin{2},varargin{3}).m;                    
                        else % relative xyz, frame
                            S.c = varargin{4}.m * anakin.vector(varargin{1},varargin{2},varargin{3}).c + varargin{4}.c; 
                            S.m = varargin{4}.m; % copy basis from given frame
                        end
                    elseif isa(varargin{4},'anakin.vector') || numel(varargin{4}) == 3 % (vector or column), ijk (vectors or columns)
                        S.c = anakin.vector(varargin{1}).c;
                        S.m = anakin.basis(varargin{2},varargin{3},varargin{4}).m;
                    else % xyz, (basis or matrix)
                        S.c = anakin.vector(varargin{1},varargin{2},varargin{3}).c;
                        S.m = anakin.basis(varargin{4}).m; 
                    end
                case 5
                    if isa(varargin{4},'anakin.basis') || numel(varargin{4}) == 9 % relative xyz, (relative basis or relative matrix), frame
                        S.c = varargin{5}.m * anakin.vector(varargin{1},varargin{2},varargin{3}).c + varargin{5}.c;
                        S.m = varargin{5}.m * anakin.basis(varargin{4}).m;
                    else % (relative vector or relative column), relative ijk, frame
                        S.c = varargin{5}.m * anakin.vector(varargin{1}).c + varargin{5}.c;
                        S.m = varargin{5}.m * anakin.basis(varargin{2},varargin{3},varargin{4}).m; 
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
