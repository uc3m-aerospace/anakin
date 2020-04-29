%{
DESCRIPTION:
returns the triangulation of a prism aligned with the Z axis and with
the geometric center at the origin. Can be used to approximate cylinders.

SYNTAX:
t = anakin.triangulations.prism(R,Lz,N)
 
AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
function t = prism(R,Lz,N)

if ~exist('R','var')
    R = 1;
end
if ~exist('Lz','var')
    Lz = 1;
end  
if ~exist('N','var')
    N = 50; % number of points to use in the contour of each basis
end  

C (4*N,3) = 0; % Allocate
C(1,:) = [2,N+1,1]; % Lower base
C(N+1,:) = [N+1,2,2*N+1]; % 1/2 of lateral
C(2*N+1,:) = [2,N+2,2*N+1]; % 1/2 of lateral
C(3*N+1,:) = [2*N+1,N+2,2*N+2]; % Upper base
for i = 2:N
    C(i,:) = [i+1,i,1]; % Lower base
    C(N+i,:) = [i,i+1,N+i]; % 1/2 of lateral
    C(2*N+i,:) = [i+1,N+i+1,N+i]; % 1/2 of lateral
    C(3*N+i,:) = [N+i,N+i+1,2*N+2]; % Upper base
end 

theta = linspace(0,2*pi*(1-1/N),N)';
x = [0;R*cos(theta);R*cos(theta);0];
y = [0;R*sin(theta);R*sin(theta);0];
z = [-Lz/2;theta*0-Lz/2;theta*0+Lz/2;+Lz/2];

t = triangulation(C,x,y,z); 
