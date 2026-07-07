function [mesh,norm,X,Y,Z] = face_creator2(center,L,nodes,angles)
%input data
C=center;
Lx=L(1);          Ly=L(2);
nx=nodes(1);      ny=nodes(2);
alfa=angles(1);   beta=angles(2);  gamma=angles(3);

%default normal versor
norm=[0,0,-1];

%mesh generation
Xinit=linspace(-Lx/2,Lx/2,nx);
Yinit=linspace(-Ly/2,Ly/2,ny);
[Xq,Yq] = meshgrid(Xinit,Yinit);
Zq = zeros(size(Yq));

%grid rotation

rot=Rot_Mat(angles);
% norm=rot*norm';
norm=norm*rot;
sz=size(Xq); 
temp=[Xq(:)-mean(Xq(:)),Yq(:)-mean(Yq(:)),Zq(:)-mean(Zq(:))]*rot;
X=reshape(temp(:,1)+mean(Xq(:)),sz)+C(1);
Y=reshape(temp(:,2)+mean(Yq(:)),sz)+C(2);
Z=reshape(temp(:,3)+mean(Zq(:)),sz)+C(3);

mesh=[X(:) Y(:) Z(:)];



% node identification vertex, edge, center
% k=1;
% for i=1:ny
%     for j=1:nx
%         mesh(k,1:3)=[X(i,j) Y(i,j) Z(i,j)];
%         if i==1 | i==ny | j==1 | j==nx
%             mesh(k,4)=2;
%         else
%             mesh(k,4)=1;
%         end
%         k=k+1;
%     end
% end
% mesh(1,4)=3;mesh(nx,4)=3;mesh(nx*ny-nx+1,4)=3;mesh(nx*ny,4)=3;
end