function [x,f,alpha_str,epsint_str,epsext_str,id_node_from,...
    id_node_to,A_link,heat,Tmin_max_st,idx]=select_trans3(x,alpha_str,epsint_str,epsext_str,...
    sat,Vf,Vf_G,v_dc,attitude,opt_st,sim_data,env,orbit,diss_mat,G_c_,...
    G_hc,f,id_node_from,id_node_to,A_link,heat,val_cases,button_handle,ROB_STRUCT)

if exist("ROB_STRUCT","var") && ~isempty(ROB_STRUCT)    
    sat = ROB_STRUCT.sat;
    sim_data = ROB_STRUCT.sim_data;
    Tmin_max_st = ROB_STRUCT.Tmin_max_st;
    T_op = ROB_STRUCT.T_op;
    n_orbit = ROB_STRUCT.n_orbit;
    dt = ROB_STRUCT.dt;
    time = ROB_STRUCT.time;
    r_dc = ROB_STRUCT.r_dc;
    rsun = ROB_STRUCT.rsun;
    sun_vers = ROB_STRUCT.sun_vers;
    a = ROB_STRUCT.a;
    Source_term_partial = ROB_STRUCT.Source_term_partial;
    if isfield(ROB_STRUCT,"Source_term_partial_cold")
        Source_term_partial_cold = ROB_STRUCT.Source_term_partial_cold;
    end

else
    Tmin_max_st=struct;
    T_op=opt_st.T_op_item;
    
    % ##### Simulation Data #####
    
    n_orbit = opt_st.n_orbit_sim;
    sim_data.n_orbit = n_orbit;
    time_step_input = opt_st.delta_t_sim;
    sim_data.step = time_step_input;
    [r_dc, ~, ~, ~, rsun, sun_vers, ~] = Orbit_analysis(env,orbit,sim_data);
    time = zeros((length(r_dc)-1)*sim_data.n_orbit,1);
    dt=r_dc(2,4)-r_dc(1,4);
    
    % ##### Environmental Heat Source Preparation #####
    
    % Environmental analysis
    [sat,a] = Environment_analysis3(sat,r_dc,v_dc,attitude,env,rsun,1);
    if isequal(val_cases,2)
        [sat,a] = Environment_analysis3(sat,r_dc,v_dc,attitude,env,rsun,2);
    end
    
    % Q source preprocessing
    [Source_term_partial,sat] = Q_source3_preprocessing(sat,env,r_dc,v_dc,attitude,rsun,a,1,Vf_G); % heat flux partial computation
    if val_cases == 2
        [Source_term_partial_cold,sat] = Q_source3_preprocessing(sat,env,r_dc,v_dc,attitude,rsun,a,2,Vf_G); % heat flux partial computation
    end
    
end

%% Solution transients

for kk=1:1:length(x(:,1)) % for every solution in the pareto front...

    if exist("button_handle","var")
        button_handle.Text = ["Transient analysis","Solution " + string(kk) + "/" + string(length(x(:,1)))];
        drawnow
    end
    % Reading optical properties input and saving them in the sat structure
    % Hypothesis: structure faces are the first 6 in sat.geom.surfaces
    for j=1:1:length(sat.geom.ext.face) % for every surface of the structure (ext)...
        % only sat.geom.surfaces is modified, sat.node.globe is not updated
        if opt_st.face_alpha_opt(j) == 1
            sat.geom.surfaces(j).prop_opt(1) = alpha_str(kk,find(opt_st.face_alpha_opt)==j);
        end
        if opt_st.face_eps_int_opt(j) == 1
            sat.geom.surfaces(j).prop_opt(2) = epsint_str(kk,find(opt_st.face_eps_int_opt)==j);
        end
        if opt_st.face_eps_ext_opt(j) == 1
            sat.geom.surfaces(j).prop_opt(3) = epsext_str(kk,find(opt_st.face_eps_ext_opt)==j);
        end
    end

    % Assigning heater power to every item specified in the optimization, for
    % the hot/default case and for the cold case if val_cases == 2
    diss_mat_new = diss_mat;
    cont = 1;
    for j = find(opt_st.item_heater_hot)
        diss_mat_new(j,2) = diss_mat_new(j,2) + heat.case1(cont,kk+1); % kk+1 because first column of heat is the id of the item (legacy data format)
        cont = cont + 1;
    end
    cont = 1;
    for j = find(opt_st.item_heater_cold)
        diss_mat_new(j,3) = diss_mat_new(j,3) + heat.case2(cont,kk+1); % kk+1 because first column of heat is the id of the item (legacy data format)
        cont = cont + 1;
    end
    % Computing nodal dissipations
    [sat]=Q_dissipation3(sat,diss_mat_new,1);
    if val_cases == 2
        [sat]=Q_dissipation3(sat,diss_mat_new,2);
    end

    % Gebhart, opt_prop.m, preparing eps_for_Ge
    eps_for_Ge = zeros(sat.node.total_node,1);
    for j=1:1:length(sat.geom.surfaces)
        node_f = sat.geom.surfaces(j).id_node;
        eps_for_Ge(node_f) = eps_for_Ge(node_f) + sat.geom.surfaces(j).prop_opt(2); % eps (int)
    end
    id_prop_eps = eps_for_Ge ~= 0; % every node that has at least 1 surface
    n_face = [sat.node.globe.n_face]';
    eps_for_Ge(id_prop_eps) = eps_for_Ge(id_prop_eps)./n_face(id_prop_eps);
    [Vf_G]=Gebhart2(Vf,eps_for_Ge,sat); % view factor correction with Gebhart method

    % Hot case heat flux
    [Source_temp_hot,sat]=Q_source3_compute(Source_term_partial,sat,r_dc,1); % heat flux final computation
    Q0_st=diag(mean(Source_temp_hot.Qa+Source_temp_hot.Qir+Source_temp_hot.Qs))+diag(sat.node.Q_diss); % average heat source for every node, performet along dimension 1 (temporal dimension) of the Q matrices
    % Cold case heat flux
    if isequal(val_cases,2)
        [Source_temp_cold,sat]=Q_source3_compute(Source_term_partial_cold,sat,r_dc,2); % heat flux final computation
        Q0_st_cold=diag(mean(Source_temp_cold.Qa_cold+Source_temp_cold.Qir_cold+Source_temp_cold.Qs_cold))+ diag(sat.node.Q_diss_cold); % average heat source for every node
    end
    
    [G0_Irr] = TMM2_only_radiative(sat,5.67e-8,Vf_G,eps_for_Ge); % computation of radiative conductor
    
    % adding conduction of thermal links
    G_c=G_c_;
    for j=1:1:opt_st.N_link
        % taxicab distance for thermal links is more conservative than euclidean distance:
        L_link = sum(abs(sat.node.globe(id_node_to(kk,j)).node - sat.node.globe(id_node_from(kk,j)).node))*10^-3; % [m]
        if L_link ~= 0
            G_link = opt_st.k_link*A_link(kk,j)*(10^(-6))/L_link;
        else
            G_link = 0;
        end
        G_c(id_node_from(kk,j),id_node_to(kk,j)) = G_c(id_node_from(kk,j),id_node_to(kk,j))+G_link;
        G_c(id_node_to(kk,j),id_node_from(kk,j)) = G_c(id_node_to(kk,j),id_node_from(kk,j,1))+G_link;
    end
    G_c = G_c-diag(sum(G_c,2)); % (?)

    %% Flag for solution validity

    % flag_1 = 1; % 0 = valid, 1 = invalid
    flag_2 = 1; % 0 = valid, 1 = invalid
    flag = 1; % 0 = valid, 1 = invalid

    %% Hot or Default Case

    % heat sources and optical properties
    Source_term.all.Qa=repmat(Source_temp_hot.Qa(1:end-1,:),n_orbit,1);
    Source_term.all.Qs=repmat(Source_temp_hot.Qs(1:end-1,:),n_orbit,1);
    Source_term.all.Qir=repmat(Source_temp_hot.Qir(1:end-1,:),n_orbit,1);
    sat.node.analysis.all.epsAf_earth=...
        repmat(sat.node.analysis.epsAf_earth(:,1:end-1)',n_orbit,1)';
    sat.node.analysis.all.epsAf_space=...
        repmat(sat.node.analysis.epsAf_space(:,1:end-1)',n_orbit,1)';
    Q_0=Source_term.all.Qs+Source_term.all.Qa+...
                Source_term.all.Qir;

    % steady-state temperature analysis
    [T_st]=Steady_comp(sat,Q0_st,G_c,G0_Irr,1);

    % transient analysis
    T_0=T_st;
    [T,time,flag_1]=transient3_pruning(sat,T_0,dt,r_dc,sim_data,time,G0_Irr,...
        G_c,G_hc,Q_0,1,T_op,opt_st);
    T=T';
    T_Celsius=T-273.15;
    Tmin_max=[min(T_Celsius,[],1); max(T_Celsius,[],1)]';
    Tmin_max_st(kk).data=Tmin_max;

    %% Cold Case

    Tmin_max_cold = [];
    T_cold_Celsius = [];
    if isequal(val_cases,2) && isequal(flag_1,0)
        
        % heat sources and optical properties
        Source_term.all.Qa_cold=repmat(Source_temp_cold.Qa_cold(1:end-1,:),n_orbit,1);
        Source_term.all.Qs_cold=repmat(Source_temp_cold.Qs_cold(1:end-1,:),n_orbit,1);
        Source_term.all.Qir_cold=repmat(Source_temp_cold.Qir_cold(1:end-1,:),n_orbit,1);
        sat.node.analysis.all.epsAf_earth_cold=...
            repmat(sat.node.analysis.epsAf_earth_cold(:,1:end-1)',n_orbit,1)';
        sat.node.analysis.all.epsAf_space_cold=...
            repmat(sat.node.analysis.epsAf_space_cold(:,1:end-1)',n_orbit,1)';
        Q_0_cold=Source_term.all.Qs_cold+Source_term.all.Qa_cold+...
            Source_term.all.Qir_cold;

        % steady-state temperature analysis
        [T_st_cold]=Steady_comp(sat,Q0_st_cold,G_c,G0_Irr,2);

        % transient analysis
        T_0=T_st_cold;
        [T_cold,time,flag_2]=transient3_pruning(sat,T_0,dt,r_dc,sim_data,time,G0_Irr,...
            G_c,G_hc,Q_0_cold,2,T_op,opt_st);
        T_cold=T_cold';
        T_cold_Celsius=T_cold-273.15;
        Tmin_max_cold=[min(T_cold_Celsius,[],1); max(T_cold_Celsius,[],1)]';
        Tmin_max_st(kk).data_cold=Tmin_max_cold;
    end

    %% Validation of the solution

    if isequal(flag_1,0) && opt_st.transient_pruning == 1
        [flag,sat] = validate_solution(sat,Tmin_max,Tmin_max_cold,val_cases,T_op,T_Celsius,T_cold_Celsius,1);
    end

    if ~exist("ROB_STRUCT","var")
        disp("Case " + string(kk) + " of " + string(length(x(:,1))) + " completd. Flag [hot, cold, final]: [" + string(flag_1) + ", " + string(flag_2) + ", " + string(flag) + "]" )
    end

    if (isequal(flag,0) || opt_st.transient_pruning == 0) % if the solution is acceptable or the pruning is disabled (not in the rob. analysis)
        if ~exist("ROB_STRUCT","var") && opt_st.save_transients == 1
            disp("Saving results of case " + string(kk))
            % load Output.mat env orbit sim_data val_cases        
            % stringa=strcat('results/trans/simcase',it,'.mat');
            % stringa=strcat('results/',erase(string(datetime("now")),[":"," "]),'trans/simcase',it,'.mat');        
            it=num2str(kk);
            stringa=strcat('results/trans/simcase',it,'_',erase(string(datetime("now")),[":"," "]),'.mat');
            sat.node.analysis.Vf=Vf;
            sat.node.analysis.Vf_G=Vf_G;
            sat.node.analysis.Orbit_data=r_dc;
            sat.node.analysis.Source_term=Source_term;
            sat.node.analysis.period=r_dc(end,4);
            sat.node.analysis.ecl_day.ecl=r_dc((r_dc(:,5)==1),1:4);
            sat.node.analysis.ecl_day.day=r_dc((r_dc(:,5)==0),1:4);
            sat.node.analysis.sun.vect=rsun;
            sat.node.analysis.sun.vers=sun_vers;
            sat.node.analysis.time_vect=time;
            save(stringa,'sat','env','orbit','sim_data','val_cases');
        end
    else
        x(kk,:)=0; % the solution is invalid
    end
end

%% Output

% alpha = alpha;
% epsint = epsint;
% epsext = epsext;
% id_node_from = id_node_from;
% id_node_to = id_node_to;
% A_link = A_link;
% heat.case1 = heat.case1;
% heat.case2 = heat.case2;

idx=find(~all(x==0,2));
x=x(idx,:);
f=f(idx,:);
if exist("alpha_str","var") && size(alpha_str,1) > 0
    alpha_str=alpha_str(idx,:);    
end
if exist("epsint_str","var") && size(epsint_str,1) > 0
    epsint_str=epsint_str(idx,:);
end
if exist("epsext_str","var") && size(epsext_str,1) > 0
    epsext_str=epsext_str(idx,:);
end
if exist("id_node_from","var") && size(id_node_from,1) > 0
    id_node_from=id_node_from(idx,:);
end
if exist("id_node_to","var") && size(id_node_to,1) > 0
    id_node_to=id_node_to(idx,:);
end
if exist("A_link","var") && size(A_link,1) > 0
    A_link=A_link(idx,:);
end
if exist("heat","var") && isfield(heat,"case1") && size(heat.case1,2) > 1
    heat.case1=heat.case1(:,[1; idx+1]);
end
if exist("heat","var") && isfield(heat,"case2") && size(heat.case2,2) > 1
    heat.case2=heat.case2(:,[1; idx+1]);
end
Tmin_max_st = Tmin_max_st(idx);

end