function [elem_cap,Con_cap] = build_sphcap(R_int,R_out,Nr_cap,Ntheta,Nt,rot,Center,z_ref,is_top,phi_block_size)

if nargin<10 || isempty(phi_block_size)
    phi_block_size=1; % vedi node_sphcap_creator.m -- 1=esatto per cella
end

h_cap = R_out; % Emisfero puro: l'unica scelta, tra le calotte a sfere
               % concentriche, che rende lo spessore della calotta
               % esattamente uguale a quello della parete (Rs_out-Rs_in =
               % R_out-R_int, verificato), oltre a essere il caso fisico
               % piu' comune per un fondo di serbatoio. Sph_Shell_Mesh.m e
               % sphcap_areas.m assumono internamente h==R_out (vedi i
               % loro commenti); h resta comunque un parametro esplicito
               % lungo tutta la catena per compatibilita' futura, se mai
               % servisse generalizzare a un cap non emisferico.

[Nodes3D,Prisms,Bricks] = Sph_Shell_Mesh(R_int,R_out,h_cap,Ntheta,Nr_cap,Nt);
Central = Tri_to_Poly(Prisms,Nt,Nr_cap+1);
total_nodes = size(Central,1)+size(Bricks,1);

% La mesh grezza ha apice a z=+R_out e bordo piatto a z=0 (calotta
% rivolta verso l'alto). Per il tappo inferiore va specchiata lungo z
% prima di posizionarla; lo specchio si fa in locale, prima della
% rotazione, poi si trasla al bordo vero (z_ref) e si applica rot+Center
% come per tutte le altre geometrie.
if ~is_top
    Nodes3D(:,3) = -Nodes3D(:,3);
end
Nodes3D(:,3) = Nodes3D(:,3) + z_ref;
Nodes3D = Nodes3D*rot + Center;

[elem_cap,Con_cap] = node_sphcap_creator(Nodes3D,Central,Bricks,R_int,R_out,h_cap,Ntheta,Nr_cap,Nt,total_nodes,phi_block_size);

% Gli ID elem.face locali (guscio interno/esterno per anello theta e
% blocco phi, celle laterali per settore) sono generati identici per il
% tappo inferiore e quello superiore da node_sphcap_creator.m, che non sa
% quale dei due sta costruendo. Qui si spostano in un blocco disgiunto per
% non collidere ne' tra loro ne' con la parete -- altrimenti surf_global.m
% mescolerebbe gusci fisicamente opposti in un'unica superficie.
%
% Lo shift deve corrispondere alla POSIZIONE reale in sat.geom.globe, non
% solo essere "disgiunto": la parete occupa sempre wall_nfaces slot
% (dimensione fissa, indipendente da phi_block_size/Ntheta -- ma dipende
% da R_int: +Nt in piu' se il cilindro e' cavo e ha anche il foro
% interno esposto, vedi cylinder_face.m/node_cyl_creator3.m), quindi il
% tappo inferiore parte subito dopo (+wall_nfaces) e quello superiore
% dopo il blocco del tappo inferiore (+wall_nfaces+block) -- vedi
% sphcap_face.m per lo stesso conteggio n_shell_slots/block usato per
% generare le normali corrispondenti.
wall_nfaces = (Nt+2) + (R_int>0)*Nt; % do_caps richiede R_int>0, quindi qui e' sempre 2*Nt+2
n_phi_blocks = ceil(Nt/phi_block_size);
n_shell_slots = 1 + (Ntheta-1)*n_phi_blocks; % apice: 1 solo slot, non n_phi_blocks -- vedi node_sphcap_creator.m
block = 2*n_shell_slots+Nt;
face_shift = wall_nfaces + is_top*block; % bottom -> +wall_nfaces, top -> +wall_nfaces+block
for idx=1:1:numel(elem_cap)
    if ~isempty(elem_cap(idx).face)
        elem_cap(idx).face = elem_cap(idx).face + face_shift;
    end
end

end
