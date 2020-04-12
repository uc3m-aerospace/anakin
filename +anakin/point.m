%{
DESCRIPTION:
point: class to model a three-dimensional geometric point. 

SYNTAX:
A0 = anakin.point();  % returns default object 
A  = anakin.point(<A|a|c>,<S1>);
where :
- <> denotes optional arguments
- | denotes alternative arguments
- () groups argument options
- A0 is the default point (origin)
- A is a point 
- a is a vector (1st-order tensor) 
- c is an array with the Cartesian coordinates
- S1 is a frame. If given, all previous input as relative to that frame
   
METHODS:  
* coordinates: returns the Cartesian coordinates of A with respect to a
  reference frame
* x: retuns a specific Cartesian coordinate of A
* displace: displace point by a vector
* pos, vel, accel: return the position, velocity and acceleration vectors
  of the point with respect to a reference frame (symbolic variables must
  be used to obtain vel and accel)   
* subs: takes values of the symbolic unknowns and returns a point with
  purely numeric components (symbolic variables must be used)   
* plot: plot a marker at the point's position

AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
classdef point 
    properties (Hidden = true, Access = protected) 
        v anakin.tensor = anakin.tensor([0;0;0]); % canonical position vector
    end
    methods % creation
        function A = point(varargin) % constructor
            switch nargin
                case 0 % no arguments
                    return;
                case 1
                    if isa(varargin{1},'anakin.point')
                        A.v = anakin.tensor(varargin{1}.v.components);
                    else
                        A.v = anakin.tensor(varargin{1});
                    end   
                case 2 % frame at the end
                    if isa(varargin{1},'anakin.point')
                        A.v = anakin.tensor(varargin{2}.origin.v.components + varargin{2}.basis.matrix * varargin{1}.v.components);
                    else
                        temp = anakin.tensor(varargin{1});
                        A.v = anakin.tensor(varargin{2}.origin.v.components + varargin{2}.basis.matrix * temp.components);
                    end                 
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in point');
            end       
        end 
        function A = set.v(A,value) % on setting v
            A.v = anakin.tensor(value);
        end
    end 
    methods (Hidden = true) % overloads
        function value = eq(A1,A2) % overload ==
            value = (A1.v == A2.v);
        end
        function value = ne(A1,A2) % overload =~
            value = ~eq(A1,A2);
        end
        function disp(A) % display
            disp('Point with coordinates:')
            disp(A.pos.components)
        end
    end
    methods % general functionality             
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
            A_.v = A.v + a;
        end    
        function pos = pos(A,S1) % Returns the position vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                pos = A.v - S1.v;
            else
                pos = A.v;
            end 
        end        
        function vel = vel(A,S1) % Returns the velocity vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                r = A.pos(S1);
                vel = r.dt(S1.basis); 
            else
                r = A.pos;
                vel = r.dt; 
            end                        
        end  
        function accel = accel(A,S1) % Returns the  acceleration vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                v = A.vel(S1);
                accel = v.dt(S1.basis);
            else
                v = A.vel;
                accel = v.dt;
            end
        end
        function A_ = subs(A,variables,values) % Particularize symbolic point
            A_ = A;   
            A_.v = A.v.subs(variables,values); 
        end          
        function h = plot(A,varargin) % plot
            c = A.v.components; 
            h = line;
            set(h,'XData',c(1),'YData',c(2),'ZData',c(3),'color','k','marker','.','linestyle','none');
            if ~isempty(varargin)
                set(h,varargin{:}); % set options stored in varargin
            end
        end
    end 
end




