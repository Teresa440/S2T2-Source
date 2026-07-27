function [elem_all,Con_all] = stitch_cyl_wall_and_caps(elem_wall,Con_wall,elem_cap_bottom,Con_cap_bottom,elem_cap_top,Con_cap_top,Nt,Nr,Nz_wall)


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
k_layer=@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;  

for j=1:1:Nr
    for i=1:1:Nt
        %% Bottom: wall h=1 <-> cap_bottom layer 
        m_wall = k_wall(i,j,1);
        m_cap  = off_b + k_layer(i,j,1);

        a_contact = elem_all(m_wall).Af(6); 
        elem_all(m_wall).Ac(6) = a_contact;
        elem_all(m_wall).Af(6) = 0;         
       
        Con_all(m_wall,m_cap) = 6;
        Con_all(m_cap,m_wall) = 3;

        %% Top: wall h=Nz_wall-1 <-> cap_top layer 
        m_wall2 = k_wall(i,j,Nz_wall-1);
        m_cap2  = off_t + k_layer(i,j,1);

        a_contact2 = elem_all(m_wall2).Af(3); 
        elem_all(m_wall2).Ac(3) = a_contact2;
        elem_all(m_wall2).Af(3) = 0;          
        
        Con_all(m_wall2,m_cap2) = 3;
        Con_all(m_cap2,m_wall2) = 6;
    end
end

end
