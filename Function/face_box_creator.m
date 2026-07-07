function [face] = face_box_creator(centre,L,n,angles)

C=centre;    
Lx=L(1);   Ly=L(2);   Lz=L(3);
nx=n(1);   ny=n(2);   nz=n(3);

r1=[1 0 0; 0 cosd(angles(1)) -sind(angles(1)); 0 sind(angles(1)) cosd(angles(1))];
r2=[cosd(angles(2)) 0 sind(angles(2)); 0 1 0; -sind(angles(2)) 0 cosd(angles(2))];
r3=[cosd(angles(3)) -sind(angles(3)) 0; sind(angles(3)) cosd(angles(3)) 0; 0 0 1];
rot=r3*r2*r1;

face=struct('ID',cell(1,6),...
            'mesh',cell(1,6),...
            'norm',cell(1,6),...
            'gridX',cell(1,6),...
            'gridY',cell(1,6),...
            'gridZ',cell(1,6));

for n=1:6
    switch n
        case 1
            angles=[90,0,0]';
            L=[Lx,Lz]';
            nodes=[nx,nz]';
            center=[0,-Ly/2,0]';
            [mesh,norm,~,~,~] = face_creator2(center,L,nodes,angles);
            face(n).ID=n;
            face(n).mesh=mesh*rot+C;
            face(n).norm=norm*rot; % removed ' from norm
            face(n).gridX=reshape(face(n).mesh(:,1),nz,[]);
            face(n).gridY=reshape(face(n).mesh(:,2),nz,[]);
            face(n).gridZ=reshape(face(n).mesh(:,3),nz,[]);

        case 2
             angles=[0,90,0]';
            L=[Lz,Ly]';
            nodes=[nz,ny]';
            center=[Lx/2,0,0]';
            [mesh,norm,~,~,~] = face_creator2(center,L,nodes,angles);
            face(n).ID=n;
            face(n).mesh=mesh*rot+C;
            face(n).norm=norm*rot; % removed ' from norm
            face(n).gridX=reshape(face(n).mesh(:,1),ny,[]);
            face(n).gridY=reshape(face(n).mesh(:,2),ny,[]);
            face(n).gridZ=reshape(face(n).mesh(:,3),ny,[]);

        case 3
            angles=[-90,0,0]';
            L=[Lx,Lz]';
            nodes=[nx,nz]';
            center=[0,Ly/2,0]';
            [mesh,norm,~,~,~] = face_creator2(center,L,nodes,angles);
            face(n).ID=n;
            face(n).mesh=mesh*rot+C;
            face(n).norm=norm*rot; % removed ' from norm
            face(n).gridX=reshape(face(n).mesh(:,1),nz,[]);
            face(n).gridY=reshape(face(n).mesh(:,2),nz,[]);
            face(n).gridZ=reshape(face(n).mesh(:,3),nz,[]);

        case 4
            angles=[0,-90,0]';
            L=[Lz,Ly]';
            nodes=[nz,ny]';
            center=[-Lx/2,0,0]';
            [mesh,norm,~,~,~] = face_creator2(center,L,nodes,angles);
            face(n).ID=n;
            face(n).mesh=mesh*rot+C;
            face(n).norm=norm*rot; % removed ' from norm
            face(n).gridX=reshape(face(n).mesh(:,1),ny,[]);
            face(n).gridY=reshape(face(n).mesh(:,2),ny,[]);
            face(n).gridZ=reshape(face(n).mesh(:,3),ny,[]);

        case 5
            angles=[0,0,0]';
            L=[Lx,Ly]';
            nodes=[nx,ny]';
            center=[0,0,-Lz/2]';
            [mesh,norm,~,~,~] = face_creator2(center,L,nodes,angles);
            face(n).ID=n;
            face(n).mesh=mesh*rot+C;
            face(n).norm=norm*rot; % removed ' from norm
            face(n).gridX=reshape(face(n).mesh(:,1),ny,[]);
            face(n).gridY=reshape(face(n).mesh(:,2),ny,[]);
            face(n).gridZ=reshape(face(n).mesh(:,3),ny,[]);

        case 6
            angles=[180,0,0]';
            L=[Lx,Ly]';
            nodes=[nx,ny]';
            center=[0,0,Lz/2]';
            [mesh,norm,~,~,~] = face_creator2(center,L,nodes,angles);
            face(n).ID=n;
            face(n).mesh=mesh*rot+C;
            face(n).norm=norm*rot; % removed ' from norm
            face(n).gridX=reshape(face(n).mesh(:,1),ny,[]);
            face(n).gridY=reshape(face(n).mesh(:,2),ny,[]);
            face(n).gridZ=reshape(face(n).mesh(:,3),ny,[]);
    end
end


end