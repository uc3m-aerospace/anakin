%{
frame: class to define reference frames

A frame can be created by passing a origin vector and the vector
basis. This is the recommended way to create the frame. Alternatively,
the components of the origin (either as an array or individually) and
the the vector basis matrix or its vectors can be given.  
By default, the code assumes the canonical reference frame S0 is 
being used in the definition. Optionally, a different frame for the
definition can be given as a last argument. Called without arguments,
frame returns the canonical reference frame S0.   
 
When a frame is given
- The origin vector is understood as relative to the given reference
  frame. 
- If the origin components are given (either as an array or
  individually), or if the vector basis matrix is given, they are also
  understood to be relative to the given reference frame and are rotated
  accordingly. 

Equality and non-equality operators have been overloaded to apply to
frames too.
 
PROPERTIES:
* O: position vector of the origin of the reference frame 
* B: vector basis of the reference frame
 
METHODS:
* rotatex,rotatey,rotatez: create simply rotated frame about one
  coordinated axis of another frame
* displace: create a frame displaced with respect to another frame
* vel, accel: return the velocity and acceleration vectors of the origin
  with respect to another reference frame (symbolic variables must be
  used)  
* omega, alpha: return the angular velocity and angular acceleration
  vectors with respect to another reference frame (symbolic variables
  must be used) 
* subs: takes values of the symbolic unknowns and returns a reference
  frame with purely numeric origin components and basis matrix (symbolic
  variables must be used)     
* plot: plots the frame with quiver

MMM20180802
%}
classdef frame
    properties 
        O anakin.vector = anakin.vector; % origin
        B anakin.basis = anakin.basis; % basis
    end
    methods
        function S = frame(varargin) % creator
            switch nargin
                case 0 % no arguments: use defaults
                    return;
                case 1 % argument is: origin || origin components || basis || basis matrix
                    if isa(varargin{1},'anakin.frame') % catch case where frame is called with a frame already
                        S = varargin{1};
                    elseif isa(varargin{1},'anakin.vector') || numel(varargin{1}) == 3 % Origin is given
                        S.O = anakin.vector(varargin{1});
                    elseif isa(varargin{1},'anakin.basis') || numel(varargin{1}) == 9  % Basis is given
                        S.B = anakin.basis(varargin{1});    
                    else
                        error('Wrong arguments in frame');
                    end
                case 2 % arguments are: O,B || origin components, frame || basis matrix, frame
                    if isa(varargin{2},'anakin.frame')
                        if isa(varargin{1},'anakin.vector') % origin, frame
                            S.O = anakin.vector(varargin{1}) + varargin{2}.O;
                            S.B = varargin{2}.B; % copy basis from given frame
                        elseif numel(varargin{1}) == 3 % Origin components, frame
                            S.O = anakin.vector(varargin{1},varargin{2}.B) + varargin{2}.O; % assume components are in the basis of given frame
                            S.B = varargin{2}.B; % copy basis from given frame
                        elseif isa(varargin{1},'anakin.basis') % Basis, frame
                            S.O = varargin{2}.O; % copy origin from given frame
                            S.B = anakin.basis(varargin{1});
                        elseif numel(varargin{1}) == 9  % Basis matrix, frame
                            S.O = varargin{2}.O; % copy origin from given frame
                            S.B = anakin.basis(varargin{1},varargin{2}.B); % assume matrix is with respect to the basis of given frame
                        else
                            error('Wrong arguments in frame');
                        end
                    else % origin, basis
                        S.O = anakin.vector(varargin{1});
                        S.B = anakin.basis(varargin{2});
                    end
                case 3 % arguments are: x,y,z || i,j,k || origin components, basis, frame || origin, basis matrix, frame || origin components, basis matrix, frame
                    if isa(varargin{3},'anakin.vector') % i,j,k
                        S.B = anakin.basis(varargin{1},varargin{2},varargin{3}); 
                    elseif isa(varargin{3},'anakin.frame') % last is frame
                        if isa(varargin{2},'anakin.basis') % origin components, basis, frame 
                            S.O = anakin.vector(varargin{1},varargin{3}.B) + varargin{3}.O;
                            S.B = varargin{2};
                        elseif isa(varargin{1},'anakin.vector') % origin, basis matrix, frame 
                            S.O = varargin{1} + varargin{3}.O;                      
                            S.B = anakin.basis(varargin{2},varargin{3}.B);
                        else % origin components, basis matrix, frame
                            S.O = anakin.vector(varargin{1},varargin{3}.B) + varargin{3}.O;
                            S.B = anakin.basis(varargin{2},varargin{3}.B);                
                        end
                    else % x,y,z
                        S.O = anakin.vector(varargin{1},varargin{2},varargin{3}); 
                    end
                case 4 % arguments are: origin (or its components) and ijk || xyz and basis (or its matrix) || xyz and frame || ijk and frame
                    if isa(varargin{4},'anakin.vector') % origin (or its components), ijk
                        S.O = anakin.vector(varargin{1});
                        S.B = anakin.basis(varargin{2},varargin{3},varargin{4});
                    elseif isa(varargin{4},'anakin.basis') % xyz, basis 
                        S.O = anakin.vector(varargin{1},varargin{2},varargin{3});
                        S.B = varargin{4};  
                    elseif isa(varargin{4},'anakin.frame') % last is frame
                        if isa(varargin{1},'anakin.vector') % ijk, frame                        
                            S.O = varargin{4}.O; % copy origin from given frame
                            S.B = anakin.basis(varargin{1},varargin{2},varargin{3});                    
                        else % xyz, frame
                            S.O = anakin.vector(varargin{1},varargin{2},varargin{3},varargin{4}.B) + varargin{4}.O;
                            S.B = varargin{4}.B; % copy basis from given frame
                        end
                    else  % xyz, basis matrix
                        S.O = anakin.vector(varargin{1},varargin{2},varargin{3});
                        S.B = anakin.basis(varargin{4});  
                    end
                case 5 % arguments are: x y z, basis, frame || x y z, basis matrix, frame || origin components, i,j,k, frame
                    if isa(varargin{4},'anakin.basis') % xyz, basis, frame
                        S.O = anakin.vector(varargin{1},varargin{2},varargin{3},varargin{5}.B) + varargin{5}.O;
                        S.B = anakin.basis(varargin{4});
                    elseif isa(varargin{4},'anakin.vector') % origin components, basis, frame
                        S.O = anakin.vector(varargin{1},varargin{5}.B) + varargin{5}.O;
                        S.B = anakin.basis(varargin{2},varargin{3},varargin{4});
                    else % xyz, basis matrix, frame
                        S.O = anakin.vector(varargin{1},varargin{2},varargin{3},varargin{5}.B) + varargin{5}.O;
                        S.B = anakin.basis(varargin{4},varargin{5}.B);                    
                    end
                case 7 % arguments are: x,y,z, i,j,k, and reference frame
                    S.O = anakin.vector(varargin{1},varargin{2},varargin{3},varargin{7}.B) + varargin{7}.O;
                    S.B = anakin.basis(varargin{4},varargin{5},varargin{6});
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in frame');
            end     
        end    
        function Sx = rotatex(S,angle) % returns rotated frame about x axis of S by angle 
            Sx = anakin.frame(S.O,S.B.rotatex(angle));
        end
        function Sy = rotatey(S,angle) % returns rotated frame about y axis of S by angle
            Sy = anakin.frame(S.O,S.B.rotatey(angle));
        end
        function Sz = rotatez(S,angle) % returns rotated frame about z axis of S by angle
            Sz = anakin.frame(S.O,S.B.rotatez(angle));
        end      
        function Sz = displace(S,OO) % returns displaced frame by vector OO
            Sz = anakin.frame(S.O + OO,S.B);
        end      
    end
    methods
        function value = eq(S1,S2) % overload ==
            value = (S1.B == S2.B && S1.O == S2.O);             
        end
        function value = ne(S1,S2) % overload ~=
            value = ~eq(S1,S2);
        end
    end
    methods
        function vO_1 = vel(S,S1) % Returns the velocity vector of O with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            rO_1 = S.O - S1.O; % relative position vector of O from S1
            vO_1 = rO_1.dt(S1.B); 
        end  
        function aO_1 = accel(S,S1) % Returns the  acceleration vector of O with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            vO_1 = S.vel(S1);
            aO_1 = vO_1.dt(S1.B);  
        end
        function omega = omega(S,S1) % Returns the angular velocity omega with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            omega = S.B.omega(S1.B); % use method from underlying bases
        end
        function alpha = alpha(S,S1) % Returns the angular acceleration vector alpha with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            alpha = S.B.alpha(S1.B); % use method from underlying bases
        end
        function S_ = subs(S,variables,values) % Particularize symbolic frame
            S_ = anakin.frame(S.O.subs(variables,values),S.B.subs(variables,values));        
        end        
    end
    methods
        function h = plot(S,varargin) % plot  
            c = S.O.components;
            m = S.B.matrix;
            hold on            
            h = quiver3([c(1),c(1),c(1)],[c(2),c(2),c(2)],[c(3),c(3),c(3)],...
            m(1,:),m(2,:),m(3,:),0,'color','k');
            hold off
            if ~isempty(varargin)
                set(h,varargin{:}); % set options stored in varargin
            end
        end
    end
end
