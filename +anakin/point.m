%{
point: class to model a geometric point. Inherits from vector.

A0 = anakin.point();  % no arguments return default object 
A  = anakin.point(A|a|c|(x,y,z),<S1>);

where :
- <> denotes optional arguments
- | denotes alternative arguments
- A0 is the default point (origin)
- A is a point 
- a is a vector
- c is an array with the three vector components
- x,y,z are the three vector components 
- S1 is a frame. If given, all previous input as relative to that frame
  
METHODS: 
* pos, vel, accel: return the velocity and acceleration vectors of the origin
  with respect to another reference frame (symbolic variables must be
  used)  
* displace: displace point by a vector
* plot: plot a dot at point's position

AUTHOR: Mario Merino <mario.merino@uc3m.es>
%}
classdef point < anakin.vector 
    methods % creation
        function A = point(varargin) % constructor
            switch nargin
                case 0 % no arguments
                    return;
                case 1 % vector or column 
                    A.c = anakin.point(varargin{1},anakin.frame).c;                         
                case 2 % (relative vector or relative column), frame
                    A.c = varargin{2}.matrix * anakin.vector(varargin{1}).c + varargin{2}.c;
                case 3 % x, y, z
                    A.c = anakin.point(varargin{1},varargin{2},varargin{3},anakin.frame).c;   
                case 4 % relative x, relative y, relative z, frame
                    A.c = varargin{4}.matrix * [varargin{1};varargin{2};varargin{3}] + varargin{4}.c;
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in point');
            end       
        end 
    end 
    methods % functionality
        function rO_1 = pos(A,S1) % Returns the position vector of the point with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            rO_1 = anakin.vector(A.c - S1.c);
        end
        function vO_1 = vel(A,S1) % Returns the velocity vector of the point with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            rO_1 = A.pos(S1);
            vO_1 = rO_1.dt(S1); 
        end  
        function aO_1 = accel(A,S1) % Returns the  acceleration vector of the point with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            vO_1 = A.vel(S1);
            aO_1 = vO_1.dt(S1);
        end
        function A_ = displace(A,a) % returns displaced point by vector a
            A_ = A;
            A_.c = A.c + a.c;
        end      
    end 
    methods % plotting
        function h = plot(A,varargin) % plot
            cc = A.c;
            hold on
            h = line(cc(1),cc(2),cc(3),'color','k','marker','.','linestyle','none');
            hold off
            if ~isempty(varargin)
                set(h,varargin{:}); % set options stored in varargin
            end
        end
    end
    methods % removed methods
        function plus(~,~)
            error('This usage of point is not permitted');
        end
        function minus(~,~)
            error('This usage of point is not permitted');
        end
        function uplus(~)
            error('This usage of point is not permitted');
        end
        function uminus(~)
            error('This usage of point is not permitted');
        end
        function times(~,~)
            error('This usage of point is not permitted');
        end
        function mtimes(~,~)
            error('This usage of point is not permitted');
        end
        function rdivide(~,~)
            error('This usage of point is not permitted');
        end
        function mrdivide(~,~)
            error('This usage of point is not permitted');
        end
        function ldivide(~,~)
            error('This usage of point is not permitted');
        end
        function mldivide(~,~)
            error('This usage of point is not permitted');
        end
        function dot(~,~)
            error('This usage of point is not permitted');
        end
        function norm(~,~)
            error('This usage of point is not permitted');
        end 
        function cross(~,~)
            error('This usage of point is not permitted');
        end
        function dir(~)
            error('This usage of point is not permitted');
        end
        function magnitude(~)
            error('This usage of point is not permitted');
        end
        function components(~,~)
            error('This usage of point is not permitted');
        end
        function x(~,~)
            error('This usage of point is not permitted');
        end
        function y(~,~)
            error('This usage of point is not permitted');
        end
        function z(~,~)
            error('This usage of point is not permitted');
        end
        function dt(~,~)
            error('This usage of point is not permitted');
        end 
        function isunitary(~)
            error('This usage of point is not permitted');
        end
        function isperpendicular(~,~)
            error('This usage of point is not permitted');
        end
        function isparallel(~,~)
            error('This usage of point is not permitted');
        end
    end
end




