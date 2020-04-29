%{
DESCRIPTION:
returns the triangulation of a sphere with the geometric center at the
origin. 

SYNTAX:
t = anakin.triangulations.sphere(R,Nphi,Ntheta)
 
AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
function t = sphere(R,Nphi,Ntheta)

if ~exist('R','var')
    R = 1;
end
if ~exist('Nphi','var')
    Nphi = 50;
end  
if ~exist('Ntheta','var')
    Ntheta = 100; % number of points to use in the contour of the basis
end  

C ((Nphi-1)*Ntheta*2,3) = 0; % Allocate 
for itheta = 1:Ntheta-1
    for iphi = 1:Nphi-1  
        C(iphi+(itheta-1)*(Nphi-1),:) = (itheta-1)*Nphi + [iphi,iphi+Nphi,iphi+1];  
        C((Nphi-1)*Ntheta+iphi+(itheta-1)*(Nphi-1),:) = (itheta-1)*Nphi + [iphi+Nphi,iphi+Nphi+1,iphi+1];  
    end
end 
for iphi = 1:Nphi-1  
    C(iphi+(Ntheta-1)*(Nphi-1),:) = [(Ntheta-1)*Nphi+iphi,iphi,(Ntheta-1)*Nphi+iphi+1];  
    C(iphi+(2*Ntheta-1)*(Nphi-1),:) = [iphi,iphi+1,(Ntheta-1)*Nphi+iphi+1];  
end

phi = linspace(-pi/2,pi/2,Nphi);
theta = linspace(0,2*pi*(1-1/Ntheta),Ntheta);
x(Nphi*Ntheta,1) = 0; % Allocate
y(Nphi*Ntheta,1) = 0; % Allocate
z(Nphi*Ntheta,1) = 0; % Allocate
for itheta = 1:Ntheta
    for iphi = 1:Nphi
        x(iphi+(itheta-1)*Nphi) = R*cos(phi(iphi)).*cos(theta(itheta));
        y(iphi+(itheta-1)*Nphi) = R*cos(phi(iphi)).*sin(theta(itheta));
        z(iphi+(itheta-1)*Nphi) = R*sin(phi(iphi));
    end
end

t = triangulation(C,x,y,z); 
