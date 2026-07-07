function [Source_term_partial,sat]=Q_source3_preprocessing(sat,env,r_dc,v_dc,attitude,rsun,a,cases,Vf_G)

% Warning, a strong hypothesis is made to speed up computation:
% Hypothesis: external nodes (structure + solar panels) are always first in sat.node.globe

Source_term_partial = struct;

N_ext_node = sum([sat.node.ext.n_node]); % structure
if sat.geom.Nsp > 0
    N_ext_node = N_ext_node + sum([sat.node.sp.n_node]); % solar panels
end
Source_term_partial.ext_nodes_id = 1:N_ext_node;

max_faces_per_node = 0;
for j=Source_term_partial.ext_nodes_id
    if length(sat.node.globe(j).face) > max_faces_per_node
        max_faces_per_node = length(sat.node.globe(j).face);
    end
end

t0=length(r_dc(:,1));
epsAf_earth_partial=zeros(t0, sat.node.total_node, max_faces_per_node);
epsAf_space_partial=zeros(t0, sat.node.total_node, max_faces_per_node);

% Environment Heat flux (hot and cold cases) initialisation
Qs_partial=zeros(t0, sat.node.total_node, max_faces_per_node); %heat flux from Sun
Qa_partial=zeros(t0, sat.node.total_node, max_faces_per_node); %heat flux from Albedo
Qir_partial=zeros(t0, sat.node.total_node, max_faces_per_node); %heat flux from Planet IR

for i=1:t0 % for every time step...
    if isempty(a)
        mat=Eci2body(sat,attitude,r_dc(i,1:3),v_dc(i,1:3),cases);
    else
        mat=Eci2body2(r_dc,a(i,:));
    end
    for j=Source_term_partial.ext_nodes_id % for every external node...
        for k=1:length(sat.node.globe(j).face)
            face_id = sat.node.globe(j).face(k);
            vertf=sat.node.globe(j).vertf(((k-1)*4+1):k*4,:); % [mm]        
            centre=mean(vertf); % by column [mm]
            switch cases
                case 1
                    F_sun=sat.node.globe(j).view_factor(k).F_sun(i,1);
                    F_p=sat.node.globe(j).view_factor(k).F_p(i,1);
                    
                case 2
                    F_sun=sat.node.globe(j).view_factor(k).F_sun_cold(i,1);
                    F_p=sat.node.globe(j).view_factor(k).F_p_cold(i,1);
                    
            end
            A_b=area_polygon2(centre,vertf(:,1),vertf(:,2),vertf(:,3))*10^-6; % [m^2] (from [mm^2]) hypothesis: vertf are in circular order
            pos=(mat'*((centre').*(10^-6)+mat*(r_dc(i,1:3)')))'; % removed round [km]

            Qir_partial(i,j,k) = A_b*env.c_ir*F_p;    % partial because it lacks the multiplication by eps_ext, which is done in another function

            % epsAf_earth_partial: partial because it lacks the multiplication by eps_ext, which is done in another function
            if strcmp(sat.node.globe(j).item,'sol')==1
                epsAf_earth_partial(i,j,k)=A_b*(F_p-sum(Vf_G(j,:))*F_p)*2; % *2 because solar panel irradiates both from top and bottom
                epsAf_space_partial(i,j,k)=A_b*(1-F_p-sum(Vf_G(j,:))*(1-F_p))*2; % *2 because solar panel irradiates both from top and bottom
            elseif strcmp(sat.node.globe(j).item,'ex')==1
                epsAf_earth_partial(i,j,k)=A_b*(F_p-(sum(Vf_G(j,:))-1)*F_p);
                epsAf_space_partial(i,j,k)=A_b*(1-F_p-(sum(Vf_G(j,:))-1)*(1-F_p));
            else
                epsAf_earth_partial(i,j,k)=A_b*F_p;
                epsAf_space_partial(i,j,k)=A_b*(1-F_p);
            end
            
            if isequal(r_dc(i,5),0) % satellite is in sunlight
                cos_solzeni=dot(pos/norm(pos),rsun/norm(rsun));
                if cos_solzeni<0
                    cos_solzeni=0;
                end
                Qa_partial(i,j,k) = A_b*env.c_s*F_p*env.c_a*cos_solzeni; % partial because it lacks the multiplication by alpha, which is done in another function
                Qs_partial(i,j,k) = A_b*env.c_s*F_sun; % partial because it lacks the multiplication by alpha, which is done in another function
            end
            
        end
    end
end

% heat flux partial
Source_term_partial.Qa_partial = Qa_partial;
Source_term_partial.Qir_partial = Qir_partial;
Source_term_partial.Qs_partial = Qs_partial;
% heat flux partial, sliced to external nodes only
Source_term_partial.Qa_partial_ext = Qa_partial(:,Source_term_partial.ext_nodes_id,:);
Source_term_partial.Qir_partial_ext = Qir_partial(:,Source_term_partial.ext_nodes_id,:);
Source_term_partial.Qs_partial_ext = Qs_partial(:,Source_term_partial.ext_nodes_id,:);
% optical properties and view factors partial
Source_term_partial.epsAf_earth_partial = epsAf_earth_partial;
Source_term_partial.epsAf_space_partial = epsAf_space_partial;
% optical properties and view factors partial, sliced to external nodes only
Source_term_partial.epsAf_earth_partial_ext = epsAf_earth_partial(:,Source_term_partial.ext_nodes_id,:);
Source_term_partial.epsAf_space_partial_ext = epsAf_space_partial(:,Source_term_partial.ext_nodes_id,:);

Source_term_partial.max_faces_per_node = max_faces_per_node;

end
