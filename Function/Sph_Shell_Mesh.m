function [Nodes3D,Prisms,Bricks] = Sph_Shell_Mesh(R_int,R_out,h,Ntheta,Nr_cap,Nt)

z_c=0;
theta_max=90; % pure hemispherical cap

Nodes_unit=Sph_Cap_Mesh(theta_max,Ntheta,Nt); 
[~,Triangles,Quads]=Sph_Cap_Mesh(1,Ntheta,Nt);
if ~isempty(Triangles)
    Triangles(:,4)=-10e10;
end
Mesh2D=[Triangles; Quads];

n=Ntheta*Nt+1;

Nodes3D=[];
Bricks=[];

for i=1:1:Nr_cap

r_rim=R_int+(i-1)*(R_out-R_int)/Nr_cap;

Nodes3D=[Nodes3D; r_rim*Nodes_unit+[0,0,z_c]];

Bricks=[Bricks; [Mesh2D+(i-1)*n Mesh2D+i*n]];

end

r_rim=R_out;
Nodes3D=[Nodes3D; r_rim*Nodes_unit+[0,0,z_c]];

A=find(sum(Bricks,2)<0);
Prisms=Bricks(A,[1 2 3 5 6 7]);
Bricks(A,:)=[];

end
