function [Nodes3D,Prisms,Bricks] = Sph_Shell_Mesh(R_int,R_out,h,Ntheta,Nr_cap,Nt)

Rs_out=(R_out^2+h^2)/(2*h);
z_c=h-Rs_out;

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
Rs=sqrt(r_rim^2+(Rs_out-h)^2);
theta_max=atan2d(r_rim,Rs_out-h);
Nodes=Sph_Cap_Mesh(theta_max,Ntheta,Nt);

Nodes3D=[Nodes3D; Rs*Nodes+[0,0,z_c]];

Bricks=[Bricks; [Mesh2D+(i-1)*n Mesh2D+i*n]];

end

r_rim=R_out;
Rs=sqrt(r_rim^2+(Rs_out-h)^2);
theta_max=atan2d(r_rim,Rs_out-h);
Nodes=Sph_Cap_Mesh(theta_max,Ntheta,Nt);
Nodes3D=[Nodes3D; Rs*Nodes+[0,0,z_c]];

A=find(sum(Bricks,2)<0);
Prisms=Bricks(A,[1 2 3 5 6 7]);
Bricks(A,:)=[];

end
