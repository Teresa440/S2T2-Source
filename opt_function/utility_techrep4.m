function [value] = utility_techrep4(x,Source_term_partial,Source_term_partial_cold,sat,Vf,r_dc,G_c,diss_mat_tot,A_node_face_avg,val_cases,opt_st)

% "x" is a matrix with size: Number_Query_Points*Number_of_Variables. 
% "dissmat" is the total dissipation matrix (3 columns for cold case,
% item_ID, dissipation power hot/default, dissipation power cold.

% global function_evaluation
% function_evaluation = function_evaluation + 1;

value = zeros(size(x,1),opt_st.N_objectives);

for i = 1:1:size(x,1)

    x_one = x(i,:);
    % deconstructing the design vector x
    alpha = x_one(opt_st.x_id_start(1):opt_st.x_id_end(1));
    eps_int = x_one(opt_st.x_id_start(2):opt_st.x_id_end(2));
    eps_ext = x_one(opt_st.x_id_start(3):opt_st.x_id_end(3));
    coord_from = x_one(opt_st.x_id_start(4):opt_st.x_id_end(4)); % size = 1 x (N_link*3)
    coord_to = x_one(opt_st.x_id_start(5):opt_st.x_id_end(5));   % size = 1 x (N_link*3)
    A_link = x_one(opt_st.x_id_start(6):opt_st.x_id_end(6));
    heater_hot_P = x_one(opt_st.x_id_start(7):opt_st.x_id_end(7));
    heater_cold_P = x_one(opt_st.x_id_start(8):opt_st.x_id_end(8));
    if opt_st.N_link > 0
        % reshaping coords
        coord_from = reshape(coord_from, 3, opt_st.N_link)'; % size = N_link x 3; columns 1, 2, 3 are x, y, z of the points, respecively
        coord_to = reshape(coord_to, 3, opt_st.N_link)';     % size = N_link x 3; columns 1, 2, 3 are x, y, z of the points, respecively
        % id nodes of the "from" and "to" list
        id_node_from_ds = dsearchn(opt_st.dsearchn_points_from,opt_st.Triang_from,coord_from);
        id_node_to_ds = dsearchn(opt_st.dsearchn_points_to,opt_st.Triang_to,coord_to);
        id_node_from = opt_st.id_node_from_list(id_node_from_ds);
        id_node_to = opt_st.id_node_to_list(id_node_to_ds);
    end
    
    %% utility_techrep
    
    % Reading optical properties input and saving them in the sat structure
    % Hypothesis: structure faces are the first 6 in sat.geom.surfaces
    cont_alpha = 1;
    cont_eps_int = 1;
    cont_eps_ext = 1;
    for j=1:1:length(sat.geom.ext.face) % for every surface of the structure (ext)...
        % only sat.geom.surfaces is modified, sat.node.globe is not updated
        if opt_st.face_alpha_opt(j) == 1
            sat.geom.surfaces(j).prop_opt(1) = alpha(cont_alpha);
            cont_alpha = cont_alpha + 1;
        end
        if opt_st.face_eps_int_opt(j) == 1
            sat.geom.surfaces(j).prop_opt(2) = eps_int(cont_eps_int);
            cont_eps_int = cont_eps_int + 1;
        end
        if opt_st.face_eps_ext_opt(j) == 1
            sat.geom.surfaces(j).prop_opt(3) = eps_ext(cont_eps_ext);
            cont_eps_ext = cont_eps_ext + 1;
        end
    end
    
    % Assigning heater power to every item specified in the optimization, for
    % the hot/default case and for the cold case if val_cases == 2
    cont_h_hot = 1;
    cont_h_cold = 1;
    for j = find(opt_st.item_heater_hot)
        diss_mat_tot(j,2) = diss_mat_tot(j,2) + heater_hot_P(cont_h_hot);
        cont_h_hot = cont_h_hot + 1;
    end
    for j = find(opt_st.item_heater_cold)
        diss_mat_tot(j,3) = diss_mat_tot(j,3) + heater_cold_P(cont_h_cold);
        cont_h_cold = cont_h_cold + 1;
    end
    % Computing nodal dissipations
    [sat]=Q_dissipation3(sat,diss_mat_tot,1);
    if val_cases == 2
        [sat]=Q_dissipation3(sat,diss_mat_tot,2);
    end
    
    % opt_prop.m, preparing eps_for_Ge
    eps_for_Ge = zeros(sat.node.total_node,1);
    for j=1:1:length(sat.geom.surfaces)
        node_f = sat.geom.surfaces(j).id_node;
        eps_for_Ge(node_f) = eps_for_Ge(node_f) + sat.geom.surfaces(j).prop_opt(2); % eps (int)
    end
    id_prop_eps = eps_for_Ge ~= 0; % every node that has at least 1 surface
    n_face = [sat.node.globe.n_face]';
    eps_for_Ge(id_prop_eps) = eps_for_Ge(id_prop_eps)./n_face(id_prop_eps);
    % Gebhart
    [Vf_G]=Gebhart2(Vf,eps_for_Ge,sat); % view factor correction with Gebhart method
    
    % Hot case heat flux
    [Source_term,sat]=Q_source3_compute(Source_term_partial,sat,r_dc,1); % heat flux final computation
    Q0_st=diag(mean(Source_term.Qa+Source_term.Qir+Source_term.Qs))+diag(sat.node.Q_diss); % average heat source for every node, performet along dimension 1 (temporal dimension) of the Q matrices
    % Cold case heat flux
    if isequal(val_cases,2)
        [Source_term_cold,sat]=Q_source3_compute(Source_term_partial_cold,sat,r_dc,2); % heat flux final computation
        Q0_st_cold=diag(mean(Source_term_cold.Qa_cold+Source_term_cold.Qir_cold+Source_term_cold.Qs_cold))+ diag(sat.node.Q_diss_cold); % average heat source for every node
    end
    
    [G0_Irr] = TMM2_only_radiative(sat,5.67e-8,Vf_G,eps_for_Ge); % computation of radiative conductor
    
    % adding conduction of thermal links
    L_link = zeros(1,opt_st.N_link);
    for j=1:1:opt_st.N_link
        % taxicab distance for thermal links is more conservative than euclidean distance:
        L_link(j) = sum(abs(sat.node.globe(id_node_to(j)).node - sat.node.globe(id_node_from(j)).node))*10^-3; % [m]
        if L_link(j) ~= 0
            G_link = opt_st.k_link*A_link(j)*(10^(-6))/L_link(j);
        else
            G_link = 0;
        end
        G_c(id_node_from(j),id_node_to(j)) = G_c(id_node_from(j),id_node_to(j))+G_link;
        G_c(id_node_to(j),id_node_from(j)) = G_c(id_node_to(j),id_node_from(j))+G_link;
    end
    
    G_c = G_c-diag(sum(G_c,2)); % (?)
    
    % steady state temperature computation
    [T_st] = Steady_comp(sat,Q0_st,G_c,G0_Irr,1);
    T_st_item = zeros(opt_st.N_item,1);    
    if isequal(val_cases,2)
        [T_st_cold] = Steady_comp(sat,Q0_st_cold,G_c,G0_Irr,2);
        T_st_item_cold = zeros(opt_st.N_item,1);
    end
    % item average temperature
    node_item=zeros(opt_st.N_item,1); % numer of nodes per i-th item
    for j=1:1:sat.node.total_node
        ID_i = sat.node.globe(j).ID_item;
        T_st_item(ID_i) = T_st_item(ID_i) + T_st(j);
        if isequal(val_cases,2)
            T_st_item_cold(ID_i) = T_st_item_cold(ID_i) + T_st_cold(j);
        end
        node_item(ID_i) = node_item(ID_i) + 1;
    end
    % structure faces average temperature
    T_st_ext_face = A_node_face_avg*T_st; % size: 6xN x Nx1 = 6x1
    if isequal(val_cases,2)
        T_st_ext_face_cold = A_node_face_avg*T_st_cold; % size: 6xN x Nx1 = 6x1
    end
    
    T_st_mean_item = T_st_item./node_item;
    T_st_mean = [T_st_ext_face; T_st_mean_item(2:end)]; % 2:end -> structure is counted as 6 faces and not 1 songle item
    if isequal(val_cases,2)
        T_st_mean_item_cold = T_st_item_cold./node_item;
        T_st_mean_cold = [T_st_ext_face_cold; T_st_mean_item_cold(2:end)]; % 2:end -> structure is counted as 6 faces and not 1 songle item
    
        it=1;
        val=abs(((opt_st.T_op_item_avg+273.15)-(T_st_mean)))';
        value(i,it)=dot(val,opt_st.k_weight); % [°C] weighted average of distance from optimal temperatures (hot case)
    
        it=it+1;
        val=abs(((opt_st.T_op_item_avg+273.15)-(T_st_mean_cold)))';
        value(i,it)=dot(val,opt_st.k_weight); % [°C] weighted average of distance from optimal temperatures (cold case)
            
        if opt_st.N_heater_hot > 0
            it=it+1;
            value(i,it)=sum(heater_hot_P); % [W] total average heater power (hot case)
        end
        if opt_st.N_heater_cold > 0
            it=it+1;
            value(i,it)=sum(heater_cold_P); % [W] total average heater power (cold case)
        end
    else
        it=1;
        val=abs(((opt_st.T_op_item_avg+273.15)-(T_st_mean)))';
        value(i,it)=dot(val,opt_st.k_weight); % [°C] weighted average of distance from optimal temperatures (default case)
    
        if opt_st.N_heater_hot > 0
            it=it+1;
            value(i,it)=sum(heater_hot_P); % [W] total average heater power (default case)
        end
    end
    if opt_st.N_link > 0
        it=it+1;
        value(i,it)=sum(A_link.*L_link); % [mm^3] total volume of links/straps material
    end
end

