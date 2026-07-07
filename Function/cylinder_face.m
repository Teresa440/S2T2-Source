function [face] = cylinder_face(Nodes3D,Nt,Nr,Nz,R,zz,an,C)
id1=size(Nodes3D,1)/Nz; % nodes per axial layer (Nt*Nr+1 solid, Nt*(Nr+1) hollow)
id2=size(Nodes3D,1);
alfa=360/Nt;
r1=[1 0 0; 0 cosd(an(1)) -sind(an(1)); 0 sind(an(1)) cosd(an(1))];
r2=[cosd(an(2)) 0 sind(an(2)); 0 1 0; -sind(an(2)) 0 cosd(an(2))];
r3=[cosd(an(3)) -sind(an(3)) 0; sind(an(3)) cosd(an(3)) 0; 0 0 1];
rot=r3*r2*r1;
face=struct( 'ID',cell(1,Nt+2),...
'norm',cell(1,Nt+2),...
'mesh',cell(1,Nt+2),...
'gridX',cell(1,Nt+2),...
'gridY',cell(1,Nt+2),...
'gridZ',cell(1,Nt+2));
% bottom face
face(1).ID=1;
face(1).norm=[0 0 1];
face(1).mesh=Nodes3D(1:id1,:);
% top face
face(2).ID=2;
face(2).norm=[0 0 -1];
face(2).mesh=Nodes3D((id2-id1+1):end,:);
% norm for lateral faces
face(3).ID=3;
face(3).norm=[1 0 0]*[cosd(-alfa/2) -sind(-alfa/2) 0; sind(-alfa/2) cosd(-alfa/2) 0; 0 0 1];
for i=1:1:Nt-1
    face(i+3).ID=i+3;
    face(i+3).norm=face(3).norm*[cosd(-alfa*i) -sind(-alfa*i) 0; sind(-alfa*i) cosd(-alfa*i) 0; 0 0 1];
end
%mesh for lateral faces
nod=zeros(Nt,2);
for i=1:1:Nt
    nod(i,1)=R*cosd(360*(i-1)/Nt);
    nod(i,2)=R*sind(360*(i-1)/Nt);
end
for i=1:1:Nt-1
for j=1:1:Nz
        c=i+2;
        face(c).mesh(j,:)=[nod(i,:) zz(j)]*rot+C;
        face(c).mesh(j+Nz,:)=[nod(i+1,:) zz(j)]*rot+C;
end
end
face(c+1).mesh=[repmat(nod(end,:),Nz,1) zz']*rot+C;
face(c+1).mesh=[face(end).mesh;[repmat(nod(1,:),Nz,1) zz']*rot+C];
for i=1:1:length(face)
    face(i).norm=face(i).norm*rot;
for j=1:1:length(face(i).mesh(:,1))
        face(i).mesh(j,:)= face(i).mesh(j,:);
end
end
end