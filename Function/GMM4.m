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
    
    for j=1:1:total_nodes
        elem(j).Af = elem(j).Af/2;
    end
    
    V_true_total = L(1)*L(2)*L(3) - (L(1)-Th(2)-Th(4))*(L(2)-Th(1)-Th(3))*(L(3)-Th(5)-Th(6));
    V_naive_total = sum([elem.V]);
    for j=1:1:total_nodes
        elem(j).V = elem(j).V*(V_true_total/V_naive_total);
    end
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
if isfield(sat.geom.cyl(i),'closed') && ~isempty(sat.geom.cyl(i).closed)
    Closed=sat.geom.cyl(i).closed;
else
    Closed=false;
end

do_caps = R_int>0 && Closed;

if do_caps
    h_cap = R; % TODO: rendere parametro utente se richiesto -- vedi
               % build_sphcap.m per la motivazione della scelta
               % emisferica (spessore calotta = spessore parete, esatto
               % solo per h_cap=R)
    if L<=2*h_cap
        error('Cylinder %d: L (%.4g) must be greater than 2*R (%.4g) for spherical end caps.',i,L,2*R);
    end
    L_wall = L-2*h_cap;
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
   
    Central=[];
    total_nodes=length(Bricks(:,1));
end

Nodes3D=Nodes3D*rot+Center;

wall_nfaces = (Nt+2) + (R_int>0)*Nt; % +Nt per il foro interno, solo se cavo --
                                      % vedi cylinder_face.m/node_cyl_creator3.m
[face] = cylinder_face(Nodes3D,Nt,Nr,Nz,R,zz,Angles,Center,R_int);
for j=1:1:wall_nfaces
    face(j).ID=face(j).ID+face_counter;
    if j <= 2 % top and bottom bases of the cylinder
        face(j).prop_opt=sat.prop.cyl(i).opt(j,:);
    else % lateral surfaces (esterne e, se presenti, interne -- stessa
         % proprieta' ottica, nessuna riga dedicata nella GUI oggi)
        face(j).prop_opt=sat.prop.cyl(i).opt(3,:);
    end
end
sat.geom.cyl(i).face=face;

[elem,Connect] = node_cyl_creator3(Nodes3D,Central,Bricks,R,L_wall,Nt,Nr,Nz,total_nodes,R_int);
n_wall = numel(elem); % per distinguere parete da calotta dopo lo stitch

if do_caps
    Ntheta_cap = Nr; % TODO: rendere parametro utente se richiesto --
                     % risoluzione angolare interna della calotta, senza
                     % analogo diretto nel cilindro (vedi Sph_Cap_Mesh.m)
    phi_block_size = 1; % quanti settori phi raggruppare sotto la stessa
                         % normale per il guscio esterno/interno (1=esatto
                         % per cella, valore piu' alto=meno superfici/
                         % costo MC ray-tracing ma normale approssimata sui
                         % settori raggruppati). Misurato su un caso di
                         % test: 1 e' sia piu' preciso (9 vs 30 gradi max)
                         % sia piu' veloce (731s vs 1036s) di 4 -- nessun
                         % motivo per usare un valore piu' alto qui.
    [elem_cb,Con_cb] = build_sphcap(R_int,R,Nr,Ntheta_cap,Nt,rot,Center,-L_wall/2,false,phi_block_size);
    [elem_ct,Con_ct] = build_sphcap(R_int,R,Nr,Ntheta_cap,Nt,rot,Center, L_wall/2,true,phi_block_size);
    [elem,Connect] = stitch_cyl_wall_and_sphcap(elem,Connect,elem_cb,Con_cb,elem_ct,Con_ct,Nt,Nr,Ntheta_cap,Nz);
end

total_nodes = numel(elem); % may include the 2 end-cap meshes
for j=1:1:total_nodes
    elem(j).ID=j+node_counter;
    elem(j).ex_in='i';
    if do_caps && j>n_wall
        elem(j).item='sphcap'; % nodi della calotta: item diverso dalla
                                % parete, altrimenti is_radial_sph &&
                                % same_sphcap in TMM2.m non scatta mai
    else
        elem(j).item='cyl';
    end
    elem(j).number=i;
    elem(j).face=elem(j).face+face_counter;
    elem(j).ID_item = ID_item_counter;
end

sat.node.cyl(i).n_node=total_nodes;
sat.node.cyl(i).elements=elem;
sat.node.globe=[sat.node.globe,elem];
sat.geom.globe=[sat.geom.globe,face];

if do_caps
    % Superfici per i due blocchi di ID riservati in build_sphcap.m
    % (bottom: +wall_nfaces, top: +wall_nfaces+block_cap). Normali
    % generate da sphcap_face.m: esatte per anello theta/blocco phi sul
    % guscio, esatte (costanti) sul bordo. Proprieta' ottiche riprese da
    % sat.prop.cyl(i).opt, stessa convenzione della parete (righe 1/2
    % basi, riga 3 laterale).
    n_phi_blocks = ceil(Nt/phi_block_size);
    n_shell_slots = 1 + (Ntheta_cap-1)*n_phi_blocks; % vedi node_sphcap_creator.m
    block_cap = 2*n_shell_slots + Nt;
    for cap_side=1:2 % 1=bottom, 2=top
        is_top = (cap_side==2);
        [cap_face] = sphcap_face(Nt,Ntheta_cap,phi_block_size,rot,Center,is_top);
        for jf=1:1:block_cap
            cap_face(jf).ID = cap_face(jf).ID + face_counter + wall_nfaces + is_top*block_cap;
            if jf<=n_shell_slots
                cap_face(jf).prop_opt = sat.prop.cyl(i).opt(cap_side,:);
            else
                cap_face(jf).prop_opt = sat.prop.cyl(i).opt(3,:);
            end
        end
        sat.geom.globe=[sat.geom.globe,cap_face];
    end
end


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

if do_caps
    face_counter=face_counter+wall_nfaces+2*block_cap; % parete + due blocchi calotta (vedi build_sphcap.m)
else
    face_counter=face_counter+wall_nfaces;
end
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