function [elem_cap,Con_cap] = build_cyl_cap(R_int,R_out,Nr,Nt,Thickness,rot,Center,z_ref,is_top)


if is_top
    zz = [z_ref, z_ref+Thickness];
else
    zz = [z_ref-Thickness, z_ref];
end

%% Annular layer (R_int..R_out), same radial grid as the wall
[Nodes_l,Tri_l,Quad_l] = Circle_Mesh(R_out,Nr,Nt,R_int);
[Nodes3D_l,~,Bricks_l] = Mesh2D_to_Mesh3D(Nodes_l,Tri_l,Quad_l,zz);
total_nodes_l = size(Bricks_l,1);
Nodes3D_l = Nodes3D_l*rot+Center;
[elem_layer,Con_layer] = node_cyl_creator3(Nodes3D_l,[],Bricks_l,R_out,Thickness,Nt,Nr,2,total_nodes_l,R_int);
n_layer = numel(elem_layer);

%% Inner core (0..R_int), solid mini-cylinder
[Nodes_c,Tri_c,Quad_c] = Circle_Mesh(R_int,Nr,Nt,0);
[Nodes3D_c,Prisms_c,Bricks_c] = Mesh2D_to_Mesh3D(Nodes_c,Tri_c,Quad_c,zz);
[Central_c] = Tri_to_Poly(Prisms_c,Nt,2);
total_nodes_c = length(Central_c(:,1))+length(Bricks_c(:,1));
Nodes3D_c = Nodes3D_c*rot+Center;
[elem_core,Con_core] = node_cyl_creator3(Nodes3D_c,Central_c,Bricks_c,R_int,Thickness,Nt,Nr,2,total_nodes_c,0);
n_core = numel(elem_core);

%% Top cap
if is_top
    for idx=1:1:n_layer
        elem_layer(idx).Ac([3 6]) = elem_layer(idx).Ac([6 3]);
        elem_layer(idx).Af([3 6]) = elem_layer(idx).Af([6 3]);
    end
    for idx=1:1:n_core
        elem_core(idx).Ac([3 6]) = elem_core(idx).Ac([6 3]);
        elem_core(idx).Af([3 6]) = elem_core(idx).Af([6 3]);
    end
end

% Un solo ID di faccia locale per l'intero disco esterno (layer anulare +
% core), stessa convenzione della base piatta della parete in
% cylinder_face.m (un'unica faccia condivisa da tutti i settori, non una
% per settore) -- prima solo il ramo is_top la impostava, quello bottom
% restava con gli ID locali grezzi di node_cyl_creator3 (sbagliato:
% avrebbe prodotto tante facce separate invece di un disco unico). Solo il
% layer: il core riceve piu' sotto un secondo ID dedicato (lato cavita').
for idx=1:1:n_layer
    elem_layer(idx).face = 1;
end
for idx=1:1:n_core
    elem_core(idx).face = 1;
end


shift_local = [0,0,(zz(2)-zz(1))/2];
shift_global = shift_local*rot;
for idx=1:1:n_layer
    elem_layer(idx).node = elem_layer(idx).node+shift_global;
end
for idx=1:1:n_core
    elem_core(idx).node = elem_core(idx).node+shift_global;
end


elem_cap = [elem_layer,elem_core];
Con_cap = zeros(n_layer+n_core);
Con_cap(1:n_layer,1:n_layer) = Con_layer;
Con_cap(n_layer+1:end,n_layer+1:end) = Con_core;


k_layer=@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;                      % R_int>0 
k_core=@(i,j,h) (Nt*(Nr-1)+1)*(h-1) + (j-2)*Nt*(j>1) + 1 + i*(j>1); % R_int==0 

idx_layer_inner = k_layer(1:Nt,1,1);      
idx_core_outer  = k_core(1:Nt,Nr,1);      

for ii=1:1:Nt
    m_layer = idx_layer_inner(ii);
    m_core  = idx_core_outer(ii)+n_layer;

    a_in = elem_cap(m_layer).Af(5); 
    elem_cap(m_layer).Ac(5) = a_in;
    elem_cap(m_layer).Af(5) = 0;

    a_out = elem_cap(m_core).Af(2); 
    elem_cap(m_core).Ac(2) = a_out;
    elem_cap(m_core).Af(2) = 0;

    
    if strcmp(elem_cap(m_core).type,'s')
        elem_cap(m_core).type = 'cq';
        elem_cap(m_core).vertf = elem_cap(m_core).vertf(1:4,:);
        elem_cap(m_core).face = elem_cap(m_core).face(1); % drop the lateral face id, keep the top/bottom one
    end

    Con_cap(m_layer,m_core) = 5; % layer looking inward, toward the core
    Con_cap(m_core,m_layer) = 2; % core looking outward, toward the layer
end

% Il core (0..R_int) tocca la CAVITA' interna del cilindro sul lato
% z_ref, non la parete (a differenza del layer, che li' tocca davvero la
% parete e resta correttamente non-radiativo) -- ma node_cyl_creator3.m,
% costruendo il tappo come un solo strato assiale (Nz=2), applica comunque
% la convenzione generica "un lato e' un confine vero (Af), l'altro si
% presume continui altrove (Ac)" -- valida per gli strati INTERNI di una
% pila vera (come la parete, Nz=9), ma qui non c'e' nessun "altrove" dal
% lato z_ref del core: e' aria (la cavita'), un confine vero quanto quello
% esterno. L'area e' gia' quella giusta (finita in Ac invece che Af per la
% stessa convenzione) -- va solo rispostata indietro, con un ID di faccia
% proprio (locale 2), distinto dal disco esterno (locale 1, condiviso con
% tutto il layer). Verificato: senza questo, il tappo emetteva verso la
% cavita' solo 0.01% (quasi nullo) mentre la cavita' emetteva verso il
% tappo il 15% -- scambio non reciproco, energeticamente sbagliato.
z_ref_idx = 3 + 3*is_top; % bottom (nessuno scambio sopra) -> 3; top (scambiato sopra) -> 6
offset_local = [0 0 Thickness];
offset_global = offset_local*rot;
for idx=1:1:n_core
    m = n_layer+idx; % elem_cap, non elem_core -- e' gia' la copia concatenata
                      % restituita dalla funzione (modificare elem_core qui
                      % non avrebbe avuto nessun effetto sul risultato finale)
    a_cav = elem_cap(m).Ac(z_ref_idx);
    elem_cap(m).Af(z_ref_idx) = a_cav;
    elem_cap(m).Ac(z_ref_idx) = 0;
    % Il lato interno non ha una propria geometria salvata da
    % node_cyl_creator3.m (elemento a strato singolo) -- lo ricavo
    % specchiando il lato esterno lungo z di uno spessore Thickness
    % (stesso profilo radiale/angolare, unica differenza reale).
    vertf_outer = elem_cap(m).vertf;
    if is_top
        vertf_inner = vertf_outer - offset_global;
    else
        vertf_inner = vertf_outer + offset_global;
    end
    elem_cap(m).vertf = [vertf_outer; vertf_inner];
    elem_cap(m).face = [1, 2]; % 1=disco esterno (condiviso col layer), 2=lato cavita' (solo core)
end

% Gli ID elem.face locali sono generati identici per il tappo inferiore e
% quello superiore -- qui si spostano in un blocco disgiunto per non
% collidere ne' tra loro ne' con la parete, stessa logica di
% build_sphcap.m. block=2 (disco esterno + lato cavita' del core). do_caps
% richiede R_int>0, quindi wall_nfaces e' sempre 2*Nt+2 (vedi
% cylinder_face.m).
wall_nfaces = 2*Nt+2;
block = 2;
face_shift = wall_nfaces + is_top*block; % bottom -> +wall_nfaces, top -> +wall_nfaces+block
for idx=1:1:numel(elem_cap)
    if ~isempty(elem_cap(idx).face)
        elem_cap(idx).face = elem_cap(idx).face + face_shift;
    end
end

end
