function [elem_cap,Con_cap] = build_sphcap(R_int,R_out,Nr_cap,Ntheta,Nt,rot,Center,z_ref,is_top)

h_cap = R_out; % TODO: rendere parametro utente se richiesto -- le formule
               % di Sph_Cap_Mesh.m/sphcap_areas.m sono gia' generali
               % (theta_max libero). h_cap=R_out (emisfero) e' l'unica
               % scelta, tra le calotte a sfere concentriche, che rende lo
               % spessore della calotta esattamente uguale a quello della
               % parete (Rs_out-Rs_in = R_out-R_int, verificato), oltre a
               % essere il caso fisico piu' comune per un fondo di
               % serbatoio.

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

[elem_cap,Con_cap] = node_sphcap_creator(Nodes3D,Central,Bricks,R_int,R_out,h_cap,Ntheta,Nr_cap,Nt,total_nodes);

% Gli ID elem.face locali (1/2 guscio esterno/interno, i+2 celle laterali,
% i=1..Nt) sono generati identici per il tappo inferiore e quello
% superiore da node_sphcap_creator.m, che non sa quale dei due sta
% costruendo. Qui si spostano in un blocco disgiunto (Nt+2 ID a testa) per
% non collidere ne' tra loro ne' con la parete (che occupa 1..Nt+2) --
% altrimenti surf_global.m mescolerebbe gusci fisicamente opposti in
% un'unica superficie. GMM4.m deve riservare 3*(Nt+2) ID invece di Nt+2
% quando do_caps e' vero, e riempire i due blocchi aggiuntivi in
% sat.geom.globe (normali/prop_opt) di conseguenza.
block = Nt+2;
face_shift = (1+is_top)*block; % bottom -> +block, top -> +2*block
for idx=1:1:numel(elem_cap)
    if ~isempty(elem_cap(idx).face)
        elem_cap(idx).face = elem_cap(idx).face + face_shift;
    end
end

% Le normali "aggregate" per questi due blocchi restano un'approssimazione
% (vedi GMM4.m, dove si costruiscono) -- una cupola non ha una normale
% unica come un disco piatto. Il calcolo dei fattori di vista sulla
% calotta e' quindi accettabile ma non fisicamente esatto finche' non si
% costruisce una generazione di facce dedicata (analoga a
% cylinder_face.m ma con normali locali per settore, non ancora fatta).

end
