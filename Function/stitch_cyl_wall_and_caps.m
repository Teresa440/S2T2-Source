function [elem_all,Con_all] = stitch_cyl_wall_and_caps(elem_wall,Con_wall,elem_cap_bottom,Con_cap_bottom,elem_cap_top,Con_cap_top,Nt,Nr,Nz_wall,Nr_cap)


n_wall  = numel(elem_wall);
n_cap_b = numel(elem_cap_bottom);
n_cap_t = numel(elem_cap_top);
n_all   = n_wall+n_cap_b+n_cap_t;

off_b = n_wall;
off_t = n_wall+n_cap_b;

elem_all = [elem_wall,elem_cap_bottom,elem_cap_top];
Con_all = zeros(n_all);
Con_all(1:n_wall,1:n_wall) = Con_wall;
Con_all(off_b+1:off_b+n_cap_b,off_b+1:off_b+n_cap_b) = Con_cap_bottom;
Con_all(off_t+1:off_t+n_cap_t,off_t+1:off_t+n_cap_t) = Con_cap_top;

k_wall =@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;
% Il layer del tappo ha risoluzione radiale FISSA (Nr_cap, oggi < Nr --
% vedi GMM4.m), indipendente da Nr della parete: piu' anelli radiali della
% parete confluiscono sullo stesso anello del tappo -- collegamento "a
% ventaglio" (mappatura j_cap sotto), non piu' 1-a-1.
k_layer=@(i,j,h) Nt*Nr_cap*(h-1) + (j-1)*Nt + i;
j_cap = @(j) min(Nr_cap, ceil(j*Nr_cap/Nr)); % anello i-esimo della parete -> anello del tappo

% L'area di contatto del tappo (Ac, facce 3/6 verso la parete) e' un
% singolo valore per elemento -- se piu' anelli della parete confluiscono
% sullo stesso anello del tappo, va divisa PRIMA del ciclo per il numero
% di anelli parete che vi confluiscono (diverso da anello a anello del
% tappo), cosi' la somma delle conduttanze in parallelo approssima la vera
% area di contatto totale invece di sovrastimarla.
if Nr_cap < Nr
    n_wall_per_cap_ring = zeros(1,Nr_cap);
    for j=1:1:Nr
        n_wall_per_cap_ring(j_cap(j)) = n_wall_per_cap_ring(j_cap(j)) + 1;
    end
    for ii=1:1:Nt
        for jc=1:1:Nr_cap
            idx_cap = k_layer(ii,jc,1);
            elem_cap_bottom(idx_cap).Ac(3) = elem_cap_bottom(idx_cap).Ac(3)/n_wall_per_cap_ring(jc);
            elem_cap_top(idx_cap).Ac(6) = elem_cap_top(idx_cap).Ac(6)/n_wall_per_cap_ring(jc);
        end
    end
    elem_all = [elem_wall,elem_cap_bottom,elem_cap_top];
end

for j=1:1:Nr
    for i=1:1:Nt
        %% Bottom: wall h=1 <-> cap_bottom layer
        m_wall = k_wall(i,j,1);
        m_cap  = off_b + k_layer(i,j_cap(j),1);

        a_contact = elem_all(m_wall).Af(6);
        elem_all(m_wall).Ac(6) = a_contact;
        elem_all(m_wall).Af(6) = 0;
        % Il disco del tappo piatto occupa lo stesso piano della faccia
        % piatta di base della parete (ID locale 2 per il tappo "Bottom",
        % verificato empiricamente -- vedi la stessa inversione trovata per
        % la calotta sferica in stitch_cyl_wall_and_sphcap.m) -- doppione
        % radiativo, va tolto da .face qui, non solo azzerato in Af/Ac.
        elem_all(m_wall).face(elem_all(m_wall).face==2) = [];
        % Nota: il lato del tappo che tocca la parete NON va toccato qui --
        % e' gia' correttamente non-radiativo per costruzione
        % (build_cyl_cap.m: Af(3)=0/Ac(3)=area reale per il tappo bottom
        % grazie alla convenzione di default di node_cyl_creator3 per
        % l'ultimo layer assiale, gia' verificato con dati reali). Un primo
        % tentativo di "correggerlo" qui aveva sovrascritto per errore
        % un'area di conduzione gia' giusta con zero -- rimosso.

        Con_all(m_wall,m_cap) = 6;
        Con_all(m_cap,m_wall) = 3;

        %% Top: wall h=Nz_wall-1 <-> cap_top layer
        m_wall2 = k_wall(i,j,Nz_wall-1);
        m_cap2  = off_t + k_layer(i,j_cap(j),1);

        a_contact2 = elem_all(m_wall2).Af(3);
        elem_all(m_wall2).Ac(3) = a_contact2;
        elem_all(m_wall2).Af(3) = 0;
        % vedi nota sopra -- stesso doppione, faccia locale 1 (base
        % superiore della parete).
        elem_all(m_wall2).face(elem_all(m_wall2).face==1) = [];

        Con_all(m_wall2,m_cap2) = 3;
        Con_all(m_cap2,m_wall2) = 6;
    end
end

end
