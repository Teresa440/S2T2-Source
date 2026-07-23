function [elem_all,Con_all] = stitch_cyl_wall_and_caps(elem_wall,Con_wall,elem_cap_bottom,Con_cap_bottom,elem_cap_top,Con_cap_top,Nt,Nr,Nz_wall)
%STITCH_CYL_WALL_AND_CAPS Combine a hollow-cylinder wall mesh with its two
%end caps (built by build_cyl_cap.m) into one elem/Con pair, adding the
%axial conduction links between the wall's h=1/h=Nz_wall-1 rings and the
%matching (same i,j) rings of each cap's own annular layer.
%
%The radial link inside each cap (layer<->core) is already wired by
%build_cyl_cap.m; this function only adds the new wall<->cap axial links.

n_wall  = numel(elem_wall);
n_cap_b = numel(elem_cap_bottom);
n_cap_t = numel(elem_cap_top);
n_all   = n_wall+n_cap_b+n_cap_t;

off_b = n_wall;          % bottom cap occupies off_b+1 .. off_b+n_cap_b
off_t = n_wall+n_cap_b;   % top cap occupies off_t+1 .. off_t+n_cap_t

elem_all = [elem_wall,elem_cap_bottom,elem_cap_top];
Con_all = zeros(n_all);
Con_all(1:n_wall,1:n_wall) = Con_wall;
Con_all(off_b+1:off_b+n_cap_b,off_b+1:off_b+n_cap_b) = Con_cap_bottom;
Con_all(off_t+1:off_t+n_cap_t,off_t+1:off_t+n_cap_t) = Con_cap_top;

k_wall =@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;  % wall's own indexing (R_int>0 branch)
k_layer=@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;  % cap's annular-layer indexing (same formula, its own local h=1)

for j=1:1:Nr
    for i=1:1:Nt
        %% Bottom: wall h=1 <-> cap_bottom layer (its only h)
        m_wall = k_wall(i,j,1);
        m_cap  = off_b + k_layer(i,j,1);

        a_contact = elem_all(m_wall).Af(6); % wall's bottom-facing area, currently exposed
        elem_all(m_wall).Ac(6) = a_contact;
        elem_all(m_wall).Af(6) = 0;         % no longer exposed: the cap covers it now
        % cap_bottom's own Ac(3) was already computed by node_cyl_creator3
        % (unused internally, since its own mesh has no h=2 layer) -- reused as-is.

        Con_all(m_wall,m_cap) = 6;
        Con_all(m_cap,m_wall) = 3;

        %% Top: wall h=Nz_wall-1 <-> cap_top layer (its only h)
        m_wall2 = k_wall(i,j,Nz_wall-1);
        m_cap2  = off_t + k_layer(i,j,1);

        a_contact2 = elem_all(m_wall2).Af(3); % wall's top-facing area, currently exposed
        elem_all(m_wall2).Ac(3) = a_contact2;
        elem_all(m_wall2).Af(3) = 0;          % no longer exposed: the cap covers it now
        % cap_top's own Ac(6) was already computed & relabeled by build_cyl_cap's
        % is_top swap (unused internally) -- reused as-is.

        Con_all(m_wall2,m_cap2) = 3;
        Con_all(m_cap2,m_wall2) = 6;
    end
end

end
