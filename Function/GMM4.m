function [sat,C] = GMM4(sat)


%% GEOMETRIC INFO - STRUCTURE (external)

Center=sat.geom.ext.center;
L=sat.geom.ext.size;
N=sat.geom.ext.nodes;
Angles=sat.geom.ext.angles;
Th=sat.geom.ext.th;

[face] = face_box_creator(Center,L,N,Angles);

for f=1:1:length(face)
    face(f).prop_opt=sat.prop.ext.opt(f,:);
end
sat.geom.ext.face=face;

[elem,total_nodes,Connect] = node_box_creator2(Center,L,N,Angles,Th);
sat.node.ext.n_node=total_nodes;
sat.node.ext.elements=elem;
for j=1:1:total_nodes
    elem(j).ID=j;
    elem(j).ex_in='e';
    elem(j).item='ex';
    elem(j).number=1;
    elem(j).ID_item = 1; % structure is always first
end
sat.node.globe=elem;
sat.geom.globe=face;
sat.node.ext.connectivity=Connect;
C=Connect;
% figure 
% spy(Connect)
node_counter=total_nodes;
face_counter=6;
ID_item_counter = 1;


%% GEOMETRIC INFO - SOLAR PANELS (external)

Nsp=sat.geom.Nsp;
if Nsp>0
for i=1:1:Nsp
ID_item_counter = ID_item_counter + 1;
Center=sat.geom.sp(i).center;
L=sat.geom.sp(i).sizes;
N=sat.geom.sp(i).nodes;
Angles=sat.geom.sp(i).angles;
Th=sat.geom.sp(i).th;

[mesh,norm,X,Y,Z] = face_creator2(Center,L,N,Angles);
sat.geom.sp(i).face.ID=i+face_counter;
sat.geom.sp(i).face.mesh=mesh;
sat.geom.sp(i).face.norm=norm; % rem '
sat.geom.sp(i).face.gridX=X;
sat.geom.sp(i).face.gridY=Y;
sat.geom.sp(i).face.gridZ=Z;
sat.geom.sp(i).face.prop_opt=sat.prop.sp(i).opt;

[elem,total_nodes,Connect] = node_face_creator3(Center,L,N,Angles,Th);

for j=1:1:total_nodes
    elem(j).ID=j+node_counter;
    elem(j).ex_in='e';
    elem(j).item='sol';
    elem(j).number=i;
    elem(j).face=i+face_counter;
    elem(j).ID_item = ID_item_counter;
end

sat.node.sp(i).n_node=total_nodes;
sat.node.sp(i).elements=elem;
sat.node.globe=[sat.node.globe,elem];
sat.geom.globe=[sat.geom.globe,sat.geom.sp(i).face];

sat.node.sp(i).connectivity=Connect;
row1=size(C,1);
row2=size(Connect,1);
col1=size(C,2);
col2=size(Connect,2);
zer1=zeros(row1,col2);
zer2=zeros(row2,col1);
C=[C,zer1;zer2,Connect];
% figure 
% spy(Connect)


node_counter=node_counter+total_nodes;
end

face_counter=face_counter+Nsp;
end

%% GEOMETRIC INFO - BOARD (internal)

Nb=sat.geom.Nb;
if Nb>0
for i=1:1:Nb
ID_item_counter = ID_item_counter + 1;
Center=sat.geom.board(i).center;
L=sat.geom.board(i).sizes;
N=sat.geom.board(i).nodes;
Angles=sat.geom.board(i).angles;
Th=sat.geom.board(i).th;

[mesh,norm,X,Y,Z] = face_creator2(Center,L,N,Angles);
sat.geom.board(i).face.ID=i+face_counter;
sat.geom.board(i).face.mesh=mesh;
sat.geom.board(i).face.norm=norm; % rem '
sat.geom.board(i).face.gridX=X;
sat.geom.board(i).face.gridY=Y;
sat.geom.board(i).face.gridZ=Z;
sat.geom.board(i).face.prop_opt=sat.prop.board(i).opt;

[elem,total_nodes,Connect] = node_face_creator3(Center,L,N,Angles,Th);

for j=1:1:total_nodes
    elem(j).ID=j+node_counter;
    elem(j).ex_in='i';
    elem(j).item='board';
    elem(j).number=i;
    elem(j).face=i+face_counter;
    elem(j).ID_item = ID_item_counter;
end

sat.node.board(i).n_node=total_nodes;
sat.node.board(i).elements=elem;
sat.node.globe=[sat.node.globe,elem];
sat.geom.globe=[sat.geom.globe,sat.geom.board(i).face];

sat.node.board(i).connectivity=Connect;
row1=size(C,1);
row2=size(Connect,1);
col1=size(C,2);
col2=size(Connect,2);
zer1=zeros(row1,col2);
zer2=zeros(row2,col1);
C=[C,zer1;zer2,Connect];
% figure 
% spy(Connect)
node_counter=node_counter+total_nodes;

end
face_counter=face_counter+Nb;
end




%% GEOMETRIC INFO - PARALLELEPIPED (internal)

NP=sat.geom.NP;
if NP>0
for i=1:1:NP
ID_item_counter = ID_item_counter + 1;
Center=sat.geom.parall(i).center;
L=sat.geom.parall(i).sizes;
N=sat.geom.parall(i).nodes;
Angles=sat.geom.parall(i).angles;
th=sat.geom.parall(i).th;


[face] = face_box_creator(Center,L,N,Angles);
for j=1:1:6
    face(j).ID=face(j).ID+face_counter;
    face(j).prop_opt=sat.prop.parall(i).opt(j,:);
end
sat.geom.parall(i).face=face;

if th>0
    % Hollow shell (uniform wall thickness on all 6 faces), same mesh
    % builder already used for the external Structure box.
    Th=repmat(th,1,6);
    [elem,total_nodes,Connect] = node_box_creator2(Center,L,N,Angles,Th);
else
    [elem,total_nodes,Connect] = node_solid_creator2(Center,L,N,Angles);
end

for j=1:1:total_nodes
    elem(j).ID=j+node_counter;
    elem(j).ex_in='i';
    elem(j).item='paral';
    elem(j).number=i;
    elem(j).face=elem(j).face+face_counter;
    elem(j).ID_item = ID_item_counter;
end

sat.node.paral(i).n_node=total_nodes;
sat.node.paral(i).elements=elem;
sat.node.globe=[sat.node.globe,elem];
sat.geom.globe=[sat.geom.globe,face];

sat.node.paral(i).connectivity=Connect;
row1=size(C,1);
row2=size(Connect,1);
col1=size(C,2);
col2=size(Connect,2);
zer1=zeros(row1,col2);
zer2=zeros(row2,col1);
C=[C,zer1;zer2,Connect];
% figure 
% spy(Connect)

node_counter=node_counter+total_nodes;
face_counter=face_counter+6;
end


end


%% GEOMETRIC INFO - CYLINDER (INTERNAL)


Nc=sat.geom.Nc;
if N>0
for i=1:1:Nc
ID_item_counter = ID_item_counter + 1;
R=sat.geom.cyl(i).R;
L=sat.geom.cyl(i).L;
Nr=sat.geom.cyl(i).Nr;
Nt=sat.geom.cyl(i).Nt;
Nz=sat.geom.cyl(i).Nz;
Angles=sat.geom.cyl(i).angles;
Center=sat.geom.cyl(i).center;
if isfield(sat.geom.cyl(i),'R_int') && ~isempty(sat.geom.cyl(i).R_int)
    R_int=sat.geom.cyl(i).R_int;
else
    R_int=0;
end
if isfield(sat.geom.cyl(i),'cap_thickness') && ~isempty(sat.geom.cyl(i).cap_thickness)
    Thickness=sat.geom.cyl(i).cap_thickness;
else
    Thickness=0;
end

do_caps = R_int>0 && Thickness>0;

if do_caps
    if L<=2*Thickness
        error('Cylinder %d: L (%.4g) must be greater than 2*cap_thickness (%.4g).',i,L,2*Thickness);
    end
    L_wall = L-2*Thickness;
else
    L_wall = L;
end

zz=linspace(-L_wall/2,L_wall/2,Nz);
r1=[1 0 0; 0 cosd(Angles(1)) -sind(Angles(1)); 0 sind(Angles(1)) cosd(Angles(1))];
r2=[cosd(Angles(2)) 0 sind(Angles(2)); 0 1 0; -sind(Angles(2)) 0 cosd(Angles(2))];
r3=[cosd(Angles(3)) -sind(Angles(3)) 0; sind(Angles(3)) cosd(Angles(3)) 0; 0 0 1];
rot=r3*r2*r1;

[Nodes, Triangles, Quads]=Circle_Mesh(R,Nr,Nt,R_int);
[Nodes3D,Prisms,Bricks] = Mesh2D_to_Mesh3D(Nodes,Triangles,Quads,zz);

if R_int==0
    [Central] = Tri_to_Poly(Prisms,Nt,Nz);
    total_nodes=length(Central(:,1))+length(Bricks(:,1));
else
    % Hollow cylinder: no central column, every ring (including the
    % inner bore ring) is meshed as a normal brick ring.
    Central=[];
    total_nodes=length(Bricks(:,1));
end

Nodes3D=Nodes3D*rot+Center;

[face] = cylinder_face(Nodes3D,Nt,Nr,Nz,R,zz,Angles,Center);
for j=1:1:Nt+2
    face(j).ID=face(j).ID+face_counter;
    if j <= 2 % top and bottom bases of the cylinder
        face(j).prop_opt=sat.prop.cyl(i).opt(j,:);
    else % lateral surfaces
        face(j).prop_opt=sat.prop.cyl(i).opt(3,:);
    end
end
sat.geom.cyl(i).face=face;

[elem,Connect] = node_cyl_creator3(Nodes3D,Central,Bricks,R,L_wall,Nt,Nr,Nz,total_nodes,R_int);

if do_caps
    [elem_cb,Con_cb] = build_cyl_cap(R_int,R,Nr,Nt,Thickness,rot,Center,-L_wall/2,false);
    [elem_ct,Con_ct] = build_cyl_cap(R_int,R,Nr,Nt,Thickness,rot,Center, L_wall/2,true);
    [elem,Connect] = stitch_cyl_wall_and_caps(elem,Connect,elem_cb,Con_cb,elem_ct,Con_ct,Nt,Nr,Nz);
end

total_nodes = numel(elem); % may include the 2 end-cap meshes
for j=1:1:total_nodes
    elem(j).ID=j+node_counter;
    elem(j).ex_in='i';
    elem(j).item='cyl';
    elem(j).number=i;
    elem(j).face=elem(j).face+face_counter;
    elem(j).ID_item = ID_item_counter;
end

sat.node.cyl(i).n_node=total_nodes;
sat.node.cyl(i).elements=elem;
sat.node.globe=[sat.node.globe,elem];
sat.geom.globe=[sat.geom.globe,face];


sat.node.cyl(i).connectivity=Connect;
row1=size(C,1);
row2=size(Connect,1);
col1=size(C,2);
col2=size(Connect,2);
zer1=zeros(row1,col2);
zer2=zeros(row2,col1);
C=[C,zer1;zer2,Connect];
% figure 
% spy(Connect)

node_counter=node_counter+total_nodes;

face_counter=face_counter+2+Nt;
end
end
disp(node_counter)

%% coord_nod

sat.node.tab_globe=struct2table(sat.node.globe);



%% global faces


B=sat.node.tab_globe.node;
for i=1:1:length(sat.geom.globe)
    %id nodes for each face
    A=sat.geom.globe(i).mesh;
    ia = ismember(A,B,"rows");
    sat.geom.globe(i).node_id=(sat.node.tab_globe.ID(ia));
    % plane coefficients
    p=sat.geom.globe(i).mesh(1,:);
    n=sat.geom.globe(i).norm;
    d = p(1)*n(1) + p(2)*n(2) + p(3)*n(3);
    sat.geom.globe(i).plane=[n(1),n(2),n(3),d];
end


sat.node.total_node = length(sat.node.globe);

% %% Diffusive nodes
% for i=1:1:length(sat.node.globe)
%     sat.node.globe(i).node_dif=mean(vertf)
% 
% 
% end