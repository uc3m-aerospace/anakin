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
    methods % overloads
        function value = eq(A1,A2) % overload ==
            value = (A1.pos0 == A2.pos0);
        end
        function value = ne(A1,A2)
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
        function rO_1 = pos(A,S1) % Returns the position vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                rO_1 = A.pos0 - S1.origin;
            else
                rO_1 = A.pos0;
            end 
        end        
        function vO_1 = vel(A,S1) % Returns the velocity vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                rO_1 = A.pos(S1);
                vO_1 = rO_1.dt(S1.basis); 
            else
                rO_1 = A.pos;
                vO_1 = rO_1.dt; 
            end                        
        end  
        function aO_1 = accel(A,S1) % Returns the  acceleration vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                vO_1 = A.vel(S1);
                aO_1 = vO_1.dt(S1.basis);
            else
                vO_1 = A.vel;
                aO_1 = vO_1.dt;
            end
        end
        function A_ = subs(A,variables,values) % Particularize symbolic frame
            A_ = A;   
            A_.pos0 = A.pos0.subs(variables,values); 
        end         
    end 
    methods % plotting
        function h = plot(A,varargin) % plot
            cc = A.pos0.components;
            hold on
            h = line(cc(1),cc(2),cc(3),'color','k','marker','.','linestyle','none');
            hold off
            if ~isempty(varargin)
                set(h,varargin{:}); % set options stored in varargin
            end
        end
    end 
end




