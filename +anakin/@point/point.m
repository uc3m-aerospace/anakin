%{
point: class to model a geometric point.

A point can be created by passing its absolute position vector in the canonical
reference frame. Alternatively, a relative position vector and a frame can be
specified. This class inherits from vector.
 
Equality and non-equality operators have been overloaded to apply to
frames too.

METHODS: 
* pos, vel, accel: return the velocity and acceleration vectors of the origin
  with respect to another reference frame (symbolic variables must be
  used)  
* displace: displace point by a vector
* plot: plot a dot at point's position

MMM20180802
%}
classdef point < anakin.vector 
    methods % creation
        function P = point(varargin) % creator
            switch nargin
                case 0 % no arguments
                    return;
                case 1 % vector or column 
                        P.c = anakin.vector(varargin{1}).c; 
                case 2 % (relative vector or relative column), frame
                    P.c = varargin{2}.matrix * anakin.vector(varargin{1}).c + varargin{2}.c;
                case 3 % x, y, z
                    P.c = [varargin{1};varargin{2};varargin{3}];
                case 4 % relative x, relative y, relative z, frame
                    P.c = varargin{4}.matrix * [varargin{1};varargin{2};varargin{3}] + varargin{4}.c;
                otherwise % other possibilities are not allowed
                    error('Wrong number of arguments in point');
            end       
        end 
    end 
    methods % functionality
        function rO_1 = pos(P,S1) % Returns the position vector of the point with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            rO_1 = anakin.vector(P.c - S1.c);
        end
        function vO_1 = vel(P,S1) % Returns the velocity vector of the point with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            rO_1 = P.pos(S1);
            vO_1 = rO_1.dt(S1); 
        end  
        function aO_1 = accel(P,S1) % Returns the  acceleration vector of the point with respect to reference frame S1
            if ~exist('S1','var') % If no S1 is given, assume the canonical reference frame
                S1 = anakin.frame;
            end
            vO_1 = P.vel(S1);
            aO_1 = vO_1.dt(S1);
        end
        function P_ = displace(P,a) % returns displaced point by vector a
            P_ = P;
            P_.c = P.c + a.c;
        end      
    end 
    methods % plotting
        function h = plot(P,varargin) % plot
            cc = P.c;
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




