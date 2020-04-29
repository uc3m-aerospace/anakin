%{
DESCRIPTION:
returns the triangulation of a box aligned with the XYZ axes and with the
geometric center at the origin. 

SYNTAX:
t = anakin.triangulations.box(Lx,Ly,Lz)
 
AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
function t = box(Lx,Ly,Lz)

if ~exist('Lx','var')
    Lx = 1;
end
if ~exist('Ly','var')
    Ly = 1;
end  
if ~exist('Lz','var')
    Lz = 1; % number of points to use in the contour of each basis
end  

C = [1,3,2;4,2,3;1,2,5;6,5,2;2,4,6;8,6,4;3,7,4;4,7,8;1,5,3;7,3,5;5,6,7;8,7,6];

x = [-1;+1;-1;+1;-1;+1;-1;+1]*Lx/2;
y = [-1;-1;+1;+1;-1;-1;+1;+1]*Ly/2;
z = [-1;-1;-1;-1;+1;+1;+1;+1]*Lz/2;

t = triangulation(C,x,y,z); 
