function [Nodes, Triangles, Quads]=Sph_Cap_Mesh(theta_max,Ntheta,Nt)

%Nodes (unit-sphere direction vectors)
%------------
c=0;
for j=1:1:Ntheta %Number of theta rings
for i=1:1:Nt %Number of phi sectors
c=Nt*(j-1)+i;
theta_j=j*theta_max/Ntheta;
phi_i=360*(i-1)/Nt;
Nodes(c,1)=sind(theta_j)*cosd(phi_i);
Nodes(c,2)=sind(theta_j)*sind(phi_i);
Nodes(c,3)=cosd(theta_j);
end
end

% Apex (theta=0), always included
Nodes(c+1,1)=0;
Nodes(c+1,2)=0;
Nodes(c+1,3)=1;

%Triangles (apex fan, ring 1 to apex)
%------------

Triangles=[];

for i=1:1:Nt-1
Triangles(i,1)=i;
Triangles(i,2)=i+1;
Triangles(i,3)=c+1;
end
Triangles(i+1,1)=i+1;
Triangles(i+1,2)=1;
Triangles(i+1,3)=c+1;

%Quads (bands between consecutive theta rings)
%------------

Quads=[];

for j=1:1:Ntheta-1
for i=1:1:Nt-1
d=Nt*(j-1)+i;
Quads(d,1)=Nt*j+i;
Quads(d,2)=Nt*j+i+1;
Quads(d,3)=Nt*(j-1)+i+1;
Quads(d,4)=Nt*(j-1)+i;
end
Quads(d+1,1)=Nt*j+i+1;
Quads(d+1,2)=Nt*j+1;
Quads(d+1,3)=Nt*(j-1)+1;
Quads(d+1,4)=Nt*(j-1)+i+1;
end

end
