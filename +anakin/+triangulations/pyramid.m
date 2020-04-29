%{
DESCRIPTION:
returns the triangulation of a pyramid aligned with the Z axis and with
the geometric center at the origin. Can be used to approximate cones.

SYNTAX:
t = anakin.triangulations.pyramid(R,Lz,N)
 
AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
function t = pyramid(R,Lz,N)

if ~exist('R','var')
    R = 1;
end
if ~exist('Lz','var')
    Lz = 1;
end  
if ~exist('N','var')
    N = 50; % number of points to use in the contour of the basis
end  

C (2*N,3) = 0; % Allocate
C(1,:) = [2,N+1,1]; % Base
C(N+1,:) = [N+1,2,N+2]; % Lateral 
for i = 2:N
    C(i,:) = [i+1,i,1]; % Base
    C(N+i,:) = [i,i+1,N+2]; % Lateral 
end 

theta = linspace(0,2*pi*(1-1/N),N)';
x = [0;R*cos(theta);0];
y = [0;R*sin(theta);0];
z = [-Lz/4;theta*0-Lz/4;+Lz/3];

t = triangulation(C,x,y,z); 
