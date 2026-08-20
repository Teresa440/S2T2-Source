function [elem_all,Con_all] = stitch_cyl_wall_and_sphcap(elem_wall,Con_wall,elem_cap_bottom,Con_cap_bottom,elem_cap_top,Con_cap_top,Nt,Nr,Ntheta,Nz_wall)


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

k_wall=@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;
k_cap =@(i,j,h) (Nt*(Ntheta-1)+1)*(h-1) + (j-2)*Nt*(j>1) + 1 + i*(j>1);

% j qui indicizza sia l'anello radiale del muro sia il layer radiale
% della calotta (Nr_cap=Nr, stessa griglia radiale per l'aggancio nodo a
% nodo -- vedi discussione Fase 4/6). Il bordo della calotta e' sempre il
% suo anello theta piu' esterno (j_theta=Ntheta), direzione 2.
for j=1:1:Nr
    for i=1:1:Nt
        %% Bottom: wall h=1 <-> cap_bottom, bordo (Ntheta) al layer j
        m_wall = k_wall(i,j,1);
        m_cap  = off_b + k_cap(i,Ntheta,j);

        a_contact = elem_all(m_wall).Af(6);
        elem_all(m_wall).Ac(6) = a_contact;
        elem_all(m_wall).Af(6) = 0;

        a_contact2 = elem_all(m_cap).Af(2);
        elem_all(m_cap).Ac(2) = a_contact2;
        elem_all(m_cap).Af(2) = 0;

        Con_all(m_wall,m_cap) = 6;
        Con_all(m_cap,m_wall) = 2;

        %% Top: wall h=Nz_wall-1 <-> cap_top, bordo (Ntheta) al layer j
        m_wall2 = k_wall(i,j,Nz_wall-1);
        m_cap2  = off_t + k_cap(i,Ntheta,j);

        a_contact3 = elem_all(m_wall2).Af(3);
        elem_all(m_wall2).Ac(3) = a_contact3;
        elem_all(m_wall2).Af(3) = 0;

        a_contact4 = elem_all(m_cap2).Af(2);
        elem_all(m_cap2).Ac(2) = a_contact4;
        elem_all(m_cap2).Af(2) = 0;

        Con_all(m_wall2,m_cap2) = 3;
        Con_all(m_cap2,m_wall2) = 2;
    end
end

end
