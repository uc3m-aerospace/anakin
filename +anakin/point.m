%{
DESCRIPTION:
point: class to model a three-dimensional geometric point. 

SYNTAX:
A = anakin.point();  % returns default object (origin of 3D space)
A = anakin.point(...,<S1>); 
where :
- <> denotes optional arguments  
- A is a point 
- ... denotes a list of one or more of the following. Later inputs can
  overwrite previous inputs:
    - point object: denotes the origin point of the reference frame
    - vector (1st-order tensor): denotes the origin point of the reference
      frame 
    - one dimensional array: denotes the origin point of the reference
      frame 
- S1 is a frame. If given, all previous inputs are relative to that frame
   
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
    properties (Hidden = true)
        r anakin.tensor = anakin.tensor([0;0;0]); % canonical position vector
    end
    methods % creation 
        function b = point(varargin) % constructor
            if isempty(varargin) % Default
                return;
            elseif length(varargin) > 1 && isa(varargin{end},'anakin.frame') % Last argin is frame
                S1 = varargin{end};
                varargin = varargin(1:end-1); 
            else % No frame is provided; use default
                S1 = anakin.frame;
            end 
            b.r = S1.r; 
            for i = 1:length(varargin) % later inputs overwrite former inputs
                temp = varargin{i};        
                if isa(temp,'anakin.point') % includes body, frame, particle as subclasses
                    b.r = anakin.tensor(S1.r.components + temp.r.components(S1)); 
                elseif isa(temp,'anakin.tensor')
                    b.r = anakin.tensor(S1.r.components + temp.components(S1));
                else % Array 
                    v_ = anakin.tensor(temp);
                    b.r = anakin.tensor(S1.r.components + v_.components(S1)); 
                end
            end 
        end 
        function A = set.r(A,value) % on setting r
            A.r = anakin.tensor(value);
            if A.r.ndims ~= 1
                error('Point position must be a vector');
            end
        end
    end 
    methods (Hidden = true) % overloads
        function value = eq(A1,A2) % overload ==
            value = (A1.r == A2.r);
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
            A_.r = A.r + a;
        end    
        function pos = pos(A,S1) % Returns the position vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                pos = A.r - S1.r;
            else
                pos = A.r;
            end 
        end        
        function vel = vel(A,S1) % Returns the velocity vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                r = A.pos(S1);
                vel = r.dt(S1); 
            else
                r = A.pos;
                vel = r.dt; 
            end                        
        end  
        function accel = accel(A,S1) % Returns the  acceleration vector of the point with respect to reference frame S1
            if exist('S1','var') % If no S1 is given, assume the canonical reference frame
                v = A.vel(S1);
                accel = v.dt(S1);
            else
                v = A.vel;
                accel = v.dt;
            end
        end
        function A_ = subs(A,variables,values) % Particularize symbolic point
            A_ = A;   
            A_.r = A.r.subs(variables,values); 
        end          
        function h = plot(A,varargin) % plot
            c = A.r.components; 
            h = line;
            set(h,'XData',c(1),'YData',c(2),'ZData',c(3),'color','r','marker','.','linestyle','none');
            for iv = 1:2:length(varargin)
                try
                    set(h,varargin{iv:iv+1});
                catch
                    % pass
                end
            end 
            set(gca,'DataAspectRatio',[1,1,1]);
        end
    end 
end




