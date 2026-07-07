function [value]=shadow3(sat,r_dc,rsun,mat,ext_surface_struct,N_ext_elem,N_ext_surf,time)

value = zeros(N_ext_elem,2); % first col = 1 -> hit with planet; second col = 1 -> hit with sun
p_start_cell = cell(N_ext_elem,2);
p_end_cell = cell(N_ext_elem,2);
hit_planet=zeros(N_ext_elem,N_ext_surf); % hit_planet(i,j) = 1 if element i is shadowed by surface j when looking to the planet
hit_sun=zeros(N_ext_elem,N_ext_surf); % hit_planet(i,j) = 1 if element i is shadowed by surface j when looking to the sun

r_dc=mat*(r_dc)';
rsun=mat*(rsun)';
norm_earth_dir=-r_dc/norm(r_dc); % planet is opposite to radial
norm_sun_dir=rsun/norm(rsun);

cont_node = 0;
for i=1:1:length(ext_surface_struct)
    for j=1:1:length(ext_surface_struct(i).elem)
        cont_node = cont_node + 1;
        % ext_surface_struct(i).elem(j); element to be checked
        p_start = ext_surface_struct(i).elem(j).center;
        % shadow from planet (intersection with back of solar panel is not computed)
        for k=1:1:N_ext_surf
            n_plane = ext_surface_struct(k).norm;
            p_plane = ext_surface_struct(k).center;
            [p_end,t] = line_plane_inters(p_start,norm_earth_dir,n_plane,p_plane);
            area_s = ext_surface_struct(k).area;
            vert_s = ext_surface_struct(k).vert;
            [flag1] = is_inside2(area_s,vert_s(:,1),vert_s(:,2),vert_s(:,3),p_end);
            if any(isnan(normal_from_points(vert_s)))
                error("ext_surface_struct(k).vert order incorrect")
            end
            flag2 = zeros(size(vert_s,1),1);
            for l=1:(size(vert_s,1)-1)
                flag2(l) = point_near_line(p_end,vert_s(l,:),vert_s(l+1,:));
            end
            flag2(end) = point_near_line(p_end,vert_s(end,:),vert_s(1,:));

            if t > 0 && flag1 == 1 && i ~= k && ~any(flag2) % ~any(flag2) = if p_end is not near any border -> NO HIT with S/C
                hit_planet(cont_node,k) = 1;
                p_start_cell{cont_node,1} = p_start;
                p_end_cell{cont_node,1} = p_end;
            end
        end
        % shadow from sun (intersection with back of solar panel is not computed)
        for k=1:1:N_ext_surf
            n_plane = ext_surface_struct(k).norm;
            p_plane = ext_surface_struct(k).center;
            [p_end,t] = line_plane_inters(p_start,norm_sun_dir,n_plane,p_plane);
            area_s = ext_surface_struct(k).area;
            vert_s = ext_surface_struct(k).vert;
            [flag1] = is_inside2(area_s,vert_s(:,1),vert_s(:,2),vert_s(:,3),p_end);

            if t > 0 && flag1 == 1 && i ~= k
                hit_sun(cont_node,k) = 1;
                p_start_cell{cont_node,2} = p_start;
                p_end_cell{cont_node,2} = p_end;
            end
        end
    end
end

value(:,1) = sum(hit_planet,2) >= 1; %hit_planet=1 if the ray from node to planet 
                                     %intersects a plane of the S/C
value(:,2) = sum(hit_sun,2) >= 1;    %hit_sun=1 if the ray from node to sun
                                     %intersects a plane of the S/C

N_node_sol = 0; % number of nodes of all solar panels                                     
for i = 7:1:(6 + sat.geom.Nsp) % first solar panel id = 7, Nsp = (N_ext_surf - 6)/2
    N_node_sol = N_node_sol + length(ext_surface_struct(i).elem);
end

% if both front and back of solar panels are in shadow (=1) then the node
% of the solar panel is in shadow, otherwise not (=0); the computation is
% solved with a binary multiplication
value((N_ext_elem-2*N_node_sol+1):(N_ext_elem-N_node_sol),:) = value((N_ext_elem-2*N_node_sol+1):(N_ext_elem-N_node_sol),:).*value((N_ext_elem-N_node_sol+1):(N_ext_elem),:);
value = value(1:(N_ext_elem-N_node_sol),:);


%% Plots

do_plot = 0;
if do_plot == 1

    figure(103)
    clf(103)
    plot_surf_elem(sat,0,103)
    for i= 1:1:size(value,1)
        if value(i,1)==1
            coord_s = p_start_cell{i,1};
            coord_e = p_end_cell{i,1};
            plot3(coord_s(1),coord_s(2),coord_s(3),'*k') % planet
            hold on
            plot3(coord_e(1),coord_e(2),coord_e(3),'ok') % planet
            hold on
        end
        if value(i,1)==1
            coord = [p_start_cell{i,1}; p_end_cell{i,1};];
            plot3(coord(:,1),coord(:,2),coord(:,3),'--k','LineWidth',0.2) % planet
            hold on
        end
    end
    title("Planet shadowed elements (t = +" + string(round(time)) + " s)")
    
    figure(104)
    clf(104)
    plot_surf_elem(sat,0,104)
    for i= 1:1:size(value,1)
        if value(i,2)==1
            coord_s = p_start_cell{i,2};
            coord_e = p_end_cell{i,2};
            plot3(coord_s(1),coord_s(2),coord_s(3),'*b') % sun
            hold on
            plot3(coord_e(1),coord_e(2),coord_e(3),'ob') % sun
            hold on
        end
        if value(i,2)==1
            coord = [p_start_cell{i,2}; p_end_cell{i,2};];
            plot3(coord(:,1),coord(:,2),coord(:,3),'--b','LineWidth',0.2) % sun
            hold on
        end
    end
    title("Sun shadowed elements (t = +" + string(round(time)) + " s)")
end

%% LOCAL FUNCTION

function flag2 = point_near_line(pt, v1, v2)
a = v1 - v2;
b = pt - v2;
line_vec = a ;%vector(start, end) # (3.5, 0, -1.5)
pnt_vec = b ;%vector(start, pnt)  # (1, 0, -1.5)
line_len = sqrt(sum(line_vec.^2)); % # 3.808
line_unitvec = line_vec/line_len; % # (0.919, 0.0, -0.394)
pnt_vec_scaled = pnt_vec/line_len; %  # (0.263, 0.0, -0.393)
t2 = dot(line_unitvec, pnt_vec_scaled); % # 0.397
if t2 < 0.0
    t2 = 0.0;
elseif t2 > 1.0
    t2 = 1.0;
end
nearest = line_vec* t2; %    # (1.388, 0.0, -0.595)
dist = sqrt(sum((nearest-pnt_vec).^2));% # 0.985
% nearest = nearest+v2;
tol = 1e-7;
if dist < tol
    flag2 = 1;
else
    flag2 = 0;
end