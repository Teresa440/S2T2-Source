function [Vf,VV] = MC_ray_tracing3(surf_for_MCRT,sat,ButtonHandle)

surface = surf_for_MCRT;
tot_nod = length(sat.node.globe(1,:));
ns = length(surface);
is_board_emitter = zeros(ns,1);
d_plane = zeros(ns,1);
n_plane = zeros(ns,3);
omega=zeros(3,1);

Vf = zeros(tot_nod,tot_nod);
rr = zeros(tot_nod,1);
VV = struct('pend',cell(1,tot_nod));

rays = 5000;

for s1 = 1:1:length(surface)
    [Rot] = rot_mat_emiss(surface(s1).norm, surface(s1).center, surface(s1).vert(1,:));
    surface(s1).Rot = Rot;
     is_board_emitter(s1) = strcmp(surface(s1).item,'board') || strcmp(surface(s1).item,'sol');
%     is_board_emitter(s1) = strcmp(surface(s1).item,'board');

    n_plane(s1,:) = surface(s1).norm;
    p_plane = surface(s1).center;
    d_plane(s1) = -dot(n_plane(s1,:),p_plane); % Plane offset parameter
end

for s1 = 1:1:length(surface)

    for i = 1:1:length(surface(s1).elem)
        id_p_start = randi(length(surface(s1).elem(i).p_start),rays,1);
        p_start = surface(s1).elem(i).p_start(id_p_start,:);
        ei = surface(s1).elem(i).ID;
        Rot = surface(s1).Rot;
        p_end_matrix = zeros(rays,ns,3); % rays * ns * 3 dimension (x y z)
        d_matrix = zeros(rays,ns); % rays * ns
        flag1_matrix = zeros(rays,ns); % rays * ns

        surface_match = surface(s1).match;
        dot_n_plane_p_start = zeros(rays,ns);
        for m = 1:1:length(surface_match)
            s2 = surface_match(m);   
            n_plane_current = n_plane(s2,:);
            n_plane_current = repmat(n_plane_current,rays,1);
            dot_n_plane_p_start(:,s2) = (dot(n_plane_current',p_start'))';

        end
        
        for r = 1:1:rays
            
            if is_board_emitter(s1)
                theta = acos(2*rand(1)-1);
                phi = 2*pi*rand(1);
                omega = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
            else
             [omega] = random_dir2(Rot,omega);
                
            end
            rr(ei) = rr(ei)+1;

            for m = 1:1:length(surface_match)
                
                s2 = surface_match(m);  

                % Line Plane Intersection: parametric line parameter t
                t = - (d_plane(s2) + dot_n_plane_p_start(r,s2))/(n_plane(s2,1)*omega(1) + n_plane(s2,2)*omega(2) + n_plane(s2,3)*omega(3));
                
                p_end_matrix(r,s2,1) = p_start(r,1) + omega(1)*t;
                p_end_matrix(r,s2,2) = p_start(r,2) + omega(2)*t;
                p_end_matrix(r,s2,3) = p_start(r,3) + omega(3)*t;

                d_matrix(r,s2) = t;

            end            
        end

        for m = 1:1:length(surface_match)
            s2 = surface_match(m);
            vert_s = surface(s2).vert;
            Rot_s2 = surface(s2).Rot;
            vert_s_rot = (Rot_s2')*(vert_s'); % 3x4 = 3x3 * 3x4
            p_end_rot_s2 = (Rot_s2')*(squeeze(p_end_matrix(:,s2,:))'); % 3 x rays = 3 x 3 * 3 x rays
            flag1_matrix(:,s2) = inpolygon(p_end_rot_s2(1,:),p_end_rot_s2(2,:),vert_s_rot(1,:),vert_s_rot(2,:));          
        end

        % eliminate rays:
        % - with intersection outside of surface (flag1_matrix == 0)
        % - from s1 to s1 (d_matrix == 0)
        % - intersecting in the wrong side (d_matrix < 0)
        id_delete = flag1_matrix == 0 | d_matrix <= 0;
        d_matrix(id_delete) = 1e20; % Inf could be used

        [min_d, s1_s2_intersec] = min(d_matrix,[],2); % ns * 1
        id_valid_rays = min_d ~= 1e20; % Rays that have at least 1 intersection with surfaces. Inf could be used (must be the same number as before)
        keep_flag = zeros(rays,ns); % rays * ns
        idx = sub2ind(size(keep_flag), find(id_valid_rays), s1_s2_intersec(id_valid_rays));
        keep_flag(idx) = 1; % rays * ns

        pe_x = squeeze(p_end_matrix(:,:,1));
        pe_y = squeeze(p_end_matrix(:,:,2));
        pe_z = squeeze(p_end_matrix(:,:,3));
        VV(ei).pend = [pe_x(keep_flag==1), pe_y(keep_flag==1), pe_z(keep_flag==1)];

        for m = 1:1:length(surface_match)
            s2_t = surface_match(m);
            Rot_s2_t = surface(s2_t).Rot;            
            id_ray = find(keep_flag(:,s2_t)); % rays * 1
            p_rec = squeeze(p_end_matrix(id_ray,s2_t,:)); %#ok<FNDSB> % n_id_ray * 3
            if ~iscolumn(p_rec)
                p_rec = p_rec'; % 3 x n_id_ray <--- n_id_ray * 3
            end
            p_rec_rot = (Rot_s2_t')*(p_rec); % 3 x n_id_ray = 3 x 3 * 3 x n_id_ray

            for j = 1:1:length(surface(s2_t).elem)
                vert_e = surface(s2_t).elem(j).vertf;
                ej = surface(s2_t).elem(j).ID;
                vert_e_rot = (Rot_s2_t')*(vert_e'); % 3x4 = 3x3 * 3x4

                flag2 = inpolygon(p_rec_rot(1,:),p_rec_rot(2,:),vert_e_rot(1,:),vert_e_rot(2,:));
                Vf(ei,ej) = Vf(ei,ej)+sum(flag2);
            end
        end

    end
    perc = s1/ns;
    if exist('ButtonHandle','var')        
        ButtonHandle.Text = "Surface " + string(s1) + "/" + string(ns);
        currentProg = min(round((size(ButtonHandle.Icon,2)-2)*(perc)),size(ButtonHandle.Icon,2)-2);    
        RGB = ButtonHandle.Icon;    
        RGB(2:end-1, 2:currentProg+1, 1) = 6/255;
        RGB(2:end-1, 2:currentProg+1, 2) = 176/255;
        RGB(2:end-1, 2:currentProg+1, 3) = 37/255;
        ButtonHandle.Icon = RGB;
        drawnow
    end
    disp(perc*100)
end

for i = 1:1:tot_nod
    if rr(i)>0
    Vf(i,:) = Vf(i,:)/rr(i);
    end
end

end

%% LOCAL FUNCTIONS

% [p_end,d] = line_plane_inters(p_start,omega,n_plane,p_plane);

% function [p_end, t] = line_plane_inters_2(p_start,omega,n_plane,p_plane)
% 
% d = -dot(n_plane,p_plane); % Plane offset parameter
% t = - (d + dot(n_plane,p_start)) / dot(n_plane,omega); % Parametric line parameter t
% p_end = p_start + omega*t; % Intersection coordinates
% 
% end






