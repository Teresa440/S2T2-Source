function [sat,a] = Environment_analysis3(sat,r_dc,v_dc,attitude,env,rsun,cases)
%% Initialization

% external surface initializaiton
N_ext_surf = (6 + sat.geom.Nsp); % 6 ext + Nsp panels surfaces
ext_surface_struct = sat.geom.surfaces(1:N_ext_surf); % 6 ext + Nsp panels surfaces
sat.geom.ext_surface = ext_surface_struct;
back_panel_surface_struct = sat.geom.surfaces(7:N_ext_surf); % back of the solar panels count as surfaces
for i=1:1:length(back_panel_surface_struct)
    back_panel_surface_struct(i).norm = -back_panel_surface_struct(i).norm;
end
% struct containing surface of satellite structure and both top and bottom of solar panels:
ext_surface_struct = [ext_surface_struct, back_panel_surface_struct]; 

N_ext_elem = 0;
for i=1:1:length(ext_surface_struct)
    N_ext_elem = N_ext_elem + length(ext_surface_struct(i).elem); % structure "ex" + solar panels "sol"
end

N_ext_node = sum([sat.node.ext.n_node]); % structure
if sat.geom.Nsp > 0
    N_ext_node = N_ext_node + sum([sat.node.sp.n_node]); % solar panels
end


% node-and-surface to element map initialization
map = zeros(length([sat.geom.ext_surface.id_node]),2);
map(:,1) = [sat.geom.ext_surface.id_node]';
elem_cont = 1;
for i=1:1:length(sat.geom.ext_surface)
    map(elem_cont:(elem_cont + length(sat.geom.ext_surface(i).id_node) - 1),2) = sat.geom.ext_surface(i).ID;
    elem_cont = elem_cont + length(sat.geom.ext_surface(i).id_node);
end
node_surf_to_elem_map = zeros(N_ext_node,N_ext_surf); % a map to convert from node and surface to element number
for i=1:1:size(map,1)
    node_surf_to_elem_map(map(i,1),map(i,2)) = i; % i-th element is the row number of map
end

% exposition matrix
norm_r=vecnorm(r_dc(:,1:3)')';
al=asin(env.Rp./norm_r);
B0=2/(7*pi)*(577/105-7*cos(al)+4/3*(cos(al)).^3-2/5*(cos(al)).^5+4/7*(cos(al)).^7);
B1=1/2*(sin(al)).^2;
B2=8/(7*pi)*(cos(al)-2*(cos(al)).^3+4*(cos(al)).^5-3*(cos(al)).^7);
B3=4/(7*pi)*(-cos(al)+40/3*(cos(al)).^3-91/3*(cos(al)).^5+18*(cos(al)).^7);
B4=8/(35*pi)*(5*cos(al)-35*(cos(al)).^3+63*(cos(al)).^5-33*(cos(al)).^7);

%% Environment View Factor Generation

%Attitude (info -> Spin; info1 -> Nadir; info2 -> Velocity) (### NO: info2 -> Nadir; info1 -> Velocity ###)
% Ext faces 1:6 -> -Y +X +Y -X -Z +Z

switch cases
    case 1
        tex_at='info1';
        tex_sp='info';
    case 2
        tex_at='info1_cold';
        tex_sp='info_cold';
end
if isfield(attitude,tex_at) % if attitude contains the field info1 or info1_cold -> Nadir-Velocity pointing
    
    for j=1:1:size(r_dc,1) % for every temporal instant...
        mat=Eci2body(sat,attitude,r_dc(j,1:3),v_dc(j,1:3),cases); %  ECI -> BODY
        inv_mat = mat'; % BODY -> ECI
        % ECI frame:
        for i = 1:1:length(sat.geom.ext_surface)
            switch cases
                case 1
                    sat.geom.ext_surface(i).att_vect(j,:)=inv_mat*(sat.geom.ext_surface(i).norm)';
                case 2
                    sat.geom.ext_surface(i).att_vect_cold(j,:)=inv_mat*(sat.geom.ext_surface(i).norm)';
            end
        end
        [shadow_mat]=shadow3(sat,r_dc(j,1:3),rsun,mat,ext_surface_struct,N_ext_elem,N_ext_surf,r_dc(j,4));
        for i=1:1:N_ext_node % Hypothesis: external nodes (structure + solar panels) are always first in sat.node.globe
            for k=1:length(sat.node.globe(i).face)
                % view factor of k-th surface of the i-th element
                [sat] = view_factor_calc(sat,i,j,k,mat,shadow_mat,r_dc,rsun,env,cases,B0,B1,B2,B3,B4,al,node_surf_to_elem_map);
            end
        end
    end
    
elseif isfield(attitude,tex_sp)  % if attitude contains the field info or info_cold -> random spin

    a0=[0 0 0];                                                            %initial angle
    rota=0.3*rand(size(r_dc,1),3);                                         %angular velocity random (maximum 5deg/s)    
    a=zeros(size(r_dc,1),3);                                               %rotation angle for each direction
    a_=rota.*r_dc(:,4);
    for j=1:1:size(r_dc,1) % for every temporal instant...
        a(j,:)=a0+a_(j,:);
        a0=a(j,:);
        mat=Eci2body2(r_dc,a(j,:)); %  ECI -> BODY
        % BODY rotations
        l1=rotz(a(j,3));
        l2=roty(a(j,2));
        l3=rotz(a(j,3));
        for i = 1:1:length(sat.geom.ext_surface)
            switch cases
                case 1
                    sat.geom.ext_surface(i).att_vect(j,:) = (l1*l2*l3*sat.geom.ext_surface(i).norm')';
                case 2
                    sat.geom.ext_surface(i).att_vect_cold(j,:) = (l1*l2*l3*sat.geom.ext_surface(i).norm')';
            end
        end
        [shadow_mat]=shadow3(sat,r_dc(j,1:3),rsun,mat,ext_surface_struct,N_ext_elem,N_ext_surf,r_dc(j,4));
        for i=1:1:N_ext_node % Hypothesis: external nodes (structure + solar panels) are always first in sat.node.globe
            for k=1:length(sat.node.globe(i).face)
                % view factor of k-th surface of the i-th element
                [sat] = view_factor_calc(sat,i,j,k,mat,shadow_mat,r_dc,rsun,env,cases,B0,B1,B2,B3,B4,al,node_surf_to_elem_map);
            end
        end
    end
end
if ~exist('a','var')
    a=[];
end
end


%% Local Functions

function [sat] = view_factor_calc(sat,i,j,k,mat,shadow_mat,r_dc,rsun,env,cases,B0,B1,B2,B3,B4,al,node_surf_to_elem_map)
    % Planet and sun view factor calculation.
    % Top and bottom surface orientation is significant in the case of
    % solar panels, wich may have the 2 opposite face of the same surfaces
    % exposed to the planet or the sun
    is_solar_panel = isequal(sat.node.globe(i).item,'sol');
    face_id = sat.node.globe(i).face(k);
    elem_id = node_surf_to_elem_map(i,face_id);
    if cases == 1
        face_att_top = sat.geom.ext_surface(face_id).att_vect(j,:);
        face_att_bottom = -sat.geom.ext_surface(face_id).att_vect(j,:); % opposite to top
    elseif cases == 2
        face_att_top = sat.geom.ext_surface(face_id).att_vect_cold(j,:);
        face_att_bottom = -sat.geom.ext_surface(face_id).att_vect_cold(j,:); % opposite to top
    end
    centre_v=sat.geom.ext_surface(face_id).center*10^-6; % [km]
    pos=(mat'*(centre_v' + round(mat*(r_dc(j,1:3)'))))'; % ECI position of the center of the face #face_id
    vers_pos=pos/norm(pos);
    % ### Planet View Factor F_p ###
    F_p_fun_1 = @(lambda) cos(lambda)/(norm(pos)/env.Rp)^2*(1-shadow_mat(elem_id,1)); % removed round(...,3)
    F_p_fun_2 = @(lambda) (...
            B0(j,1)+...
            B1(j,1)*cos(lambda)+...
            B2(j,1)*cos(lambda)^2+...
            B3(j,1)*cos(lambda)^4+...
            B4(j,1)*cos(lambda)^6)*(1-shadow_mat(elem_id,1)); % removed round(...,3)
    % Top surface
    lambda_top=acos(dot(-vers_pos,face_att_top));
    if lambda_top+al(j,1)<=pi/2
        F_p_top =  F_p_fun_1(lambda_top);
    else
        F_p_top =  F_p_fun_2(lambda_top);
    end
    % Bottom surface (solar panels only)
    if is_solar_panel == 1
        lambda_bottom=acos(dot(-vers_pos,face_att_bottom));
        if lambda_bottom+al(j,1)<=pi/2
            F_p_bottom =  F_p_fun_1(lambda_bottom);
        else
            F_p_bottom =  F_p_fun_2(lambda_bottom);
        end
    end
    % ### Sun View Factor F_sun ###
    % Top surface
    if r_dc(j,5) == 0 % r_dc(j,5) contains 1 if satellite is in eclipse, 0 otherwise
        F_sun_v=dot(face_att_top,(rsun)/norm(rsun))*(1-shadow_mat(elem_id,2));
        if F_sun_v>0
            F_sun_top=F_sun_v;
        else
            F_sun_top=0;
        end
    else % satellite is in eclipse
        F_sun_top=0;
    end
    % Bottom surface (solar panels only)
    if is_solar_panel == 1
        if r_dc(j,5) == 0 % r_dc(j,5) contains 1 if satellite is in eclipse, 0 otherwise
            F_sun_v=dot(face_att_bottom,(rsun)/norm(rsun))*(1-shadow_mat(elem_id,2));
            if F_sun_v>0
                F_sun_bottom=F_sun_v;
            else
                F_sun_bottom=0;
            end
        else % satellite is in eclipse
            F_sun_bottom=0;
        end
    end
    % ### saving view factor in the sat struct ###
    if is_solar_panel == 1
        % in case of a solar panel heat is absorbed from sun and planet
        % either from a face or the opposite (the one with max F)
        F_sun = max(F_sun_top, F_sun_bottom);
        F_p = max(F_p_top, F_p_bottom);
    else
        F_sun = F_sun_top;
        F_p = F_p_top;
    end
    if cases == 1
        sat.node.globe(i).view_factor(k).F_sun(j,1) = F_sun;
        sat.node.globe(i).view_factor(k).F_p(j,1) = F_p;
    elseif cases == 2
        sat.node.globe(i).view_factor(k).F_sun_cold(j,1) = F_sun;
        sat.node.globe(i).view_factor(k).F_p_cold(j,1) = F_p;
    end
end