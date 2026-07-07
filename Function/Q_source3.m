function [Source_term,sat]=Q_source3(Source_term,sat,env,r_dc,v_dc,attitude,rsun,a,cases,Vf_G)
t0=length(r_dc(:,1));

N_ext_node = sum([sat.node.ext.n_node]); % structure
if sat.geom.Nsp > 0
    N_ext_node = N_ext_node + sum([sat.node.sp.n_node]); % solar panels
end

epsAf_earth=zeros(sat.node.total_node,t0);
epsAf_space=zeros(sat.node.total_node,t0);

%Environment Heat flux (hot and cold cases) initialisation
Qs=zeros(t0,sat.node.total_node);                                          %heat flux from Sun
Qa=zeros(t0,sat.node.total_node);                                          %heat flux from Albedo
Qir=zeros(t0,sat.node.total_node);                                         %heat flux from planet

%heat flux from IR earth

for i=1:t0
    %heat flux incoming from the space
    if isempty(a)
        mat=Eci2body(sat,attitude,r_dc(i,1:3),v_dc(i,1:3),cases);
    else
        mat=Eci2body2(r_dc,a(i,:));
    end
    
    for j=1:N_ext_node % Hypothesis: external nodes (structure + solar panels) are always first in sat.node.globe
        Qs(i,j)=0;
        Qa(i,j)=0;
        Qir(i,j)=0;
        for k=1:length(sat.node.globe(j).face)
            face_id=sat.node.globe(j).face(k);            
            vertf=sat.node.globe(j).vertf(((k-1)*4+1):k*4,:); % [mm]        
            centre=mean(vertf); % by column [mm]
            % optical prop
            alpha=sat.geom.surfaces(face_id).prop_opt(1);
            if isequal(sat.node.globe(j).item,'ex') % satellite structure
                eps_ext=sat.geom.surfaces(face_id).prop_opt(3);
            elseif isequal(sat.node.globe(j).item,'sol') % satellite solar panels
                eps_ext=sat.geom.surfaces(face_id).prop_opt(2);
            end
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
            Qir(i,j)=Qir(i,j)+eps_ext*A_b*env.c_ir*F_p;

            if strcmp(sat.node.globe(j).item,'sol')==1
                epsAf_earth(j,i)=epsAf_earth(j,i)+eps_ext*A_b*(F_p-sum(Vf_G(j,:))*F_p)*2;
                epsAf_space(j,i)=epsAf_space(j,i)+eps_ext*A_b*(1-F_p-sum(Vf_G(j,:))*(1-F_p))*2;
            elseif strcmp(sat.node.globe(j).item,'ex')==1
                epsAf_earth(j,i)=epsAf_earth(j,i)+eps_ext*A_b*(F_p-(sum(Vf_G(j,:))-1)*F_p);
                epsAf_space(j,i)=epsAf_space(j,i)+eps_ext*A_b*(1-F_p-(sum(Vf_G(j,:))-1)*(1-F_p));
            else
                epsAf_earth(j,i)=epsAf_earth(j,i)+eps_ext*A_b*F_p;
                epsAf_space(j,i)=epsAf_space(j,i)+eps_ext*A_b*(1-F_p);
            end

            if isequal(r_dc(i,5),1) % satellite is in eclipse
                Qa(i,j)=Qa(i,j)+0;
                Qs(i,j)=Qs(i,j)+0;
            else
                cos_solzeni=dot(pos/norm(pos),rsun/norm(rsun));
                if cos_solzeni<0
                    cos_solzeni=0;
                end
                Qa(i,j)=Qa(i,j)+alpha*A_b*env.c_s*F_p*env.c_a*cos_solzeni;
                Qs(i,j)=Qs(i,j)+alpha*A_b*env.c_s*F_sun;
            end
        end
    end
    
    switch cases
        case 1
            Source_term.Qa=Qa;
            Source_term.Qir=Qir;
            Source_term.Qs=Qs;
            sat.node.analysis.epsAf_earth=epsAf_earth;
            sat.node.analysis.epsAf_space=epsAf_space;
        case 2
            Source_term.Qa_cold=Qa;
            Source_term.Qir_cold=Qir;
            Source_term.Qs_cold=Qs;
            sat.node.analysis.epsAf_earth_cold=epsAf_earth;
            sat.node.analysis.epsAf_space_cold=epsAf_space;
    end
end