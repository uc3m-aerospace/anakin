%{
DESCRIPTION:
point: class to model a geometric point. Inherits from vector.

SYNTAX:
A0 = anakin.point();  % returns default object 
A  = anakin.point(A|a|c,<S1>);
where :
- <> denotes optional arguments
- | denotes alternative arguments
- A0 is the default point (origin)
- A is a point 
- a is a vector
- c is an array with the coordinates of the point
- S1 is a frame. If given, all previous input as relative to that frame
   
METHODS: 
* spacedim: returns dimensionality of space
* pos, vel, accel: return the position, velocity and acceleration vectors
  of the point with respect to another reference frame (symbolic variables
  must be used in vel and accel)  
* displace: displace point by a vector
* plot: plot a dot at point's position

AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef point 
    properties (Hidden = true, Access = protected) 
        pos0 anakin.tensor = anakin.tensor([0;0;0]); 
    end
    methods % creation
        function A = point(varargin) % constructor
            switch nargin
                case 0 % no arguments
                    return;
                case 1 % vector or column 
                    A.pos0 = anakin.point(varargin{1},anakin.frame).pos0;
                case 2 % (relative vector or relative column), frame
                    A.pos0 = varargin{2}.basis.matrix * anakin.tensor(varargin{1}).components + varargin{2}.origin.coordinates; 
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in point');
            end       
        end 
    end 
    methods (Hidden = true) % overloads
        function value = eq(A1,A2) % overload ==
            value = (A1.pos0 == A2.pos0);
        end
        function value = ne(A1,A2) % overload =~
            value = ~eq(A1,A2);
        end
        function disp(A) % display
            disp('Point with canonical position vector components:')
            disp(A.pos0.components)
        end
    end
    methods % general functionality 
        function spacedim = spacedim(A) % number of dimensions of space
            spacedim = A.pos0.spacedim; 
        end            
        function coordinates = coordinates(A,S1) % returns the coordinates of A with respect to S1
            if exist('S1','var')
                coordinates = A.pos(S1).components(S1.basis);
            else
                coordinates = A.pos.components;
            end           
        end
        function x = x(A,i,S1) % returns single coordinate x(i) with respect to S1
            if exist('S1','var')
                coordinates = A.coordinates(S1);
            else
                coordinates = A.coordinates;
            end            
            x = coordinates(i);
        end
        function A_ = displace(A,a) % returns displaced point by vector a
            A_ = A;
            A_.pos0 = A.pos0 + a;
        end    
        function r = pos(A,S1) % Returns the position vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                r = A.pos0 - S1.origin;
            else
                r = A.pos0;
            end 
        end        
        function v = vel(A,S1) % Returns the velocity vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                r = A.pos(S1);
                v = r.dt(S1.basis); 
            else
                r = A.pos;
                v = r.dt; 
            end                        
        end  
        function a = accel(A,S1) % Returns the  acceleration vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                v = A.vel(S1);
                a = v.dt(S1.basis);
            else
                v = A.vel;
                a = v.dt;
            end
        end
        function A_ = subs(A,variables,values) % Particularize symbolic frame
            A_ = A;   
            A_.pos0 = A.pos0.subs(variables,values); 
        end         
    end 
    methods % plotting
        function h = plot(A,varargin) % plot
            if A.spacedim ~= 3
                error('This functionality is only available for points in 3D space');
            end
            c = A.pos0.components; 
            h = line;
            set(h,'XData',c(1),'YData',c(2),'ZData',c(3),'color','k','marker','.','linestyle','none');
            if ~isempty(varargin)
                set(h,varargin{:}); % set options stored in varargin
            end
        end
    end 
end




