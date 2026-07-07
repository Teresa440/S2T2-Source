function [opt_st, max_generation, sat, lb, ub] = setup_var3(...
    opt_st, max_fun_eval, population_size, sat,...
    alpha_lb, eps_int_lb, eps_ext_lb, A_link_lb, heater_hot_lb, heater_cold_lb,...
    alpha_ub, eps_int_ub, eps_ext_ub, A_link_ub, heater_hot_ub, heater_cold_ub, val_cases)

%% Bounds and Optimization Struct opt_st

max_generation = ceil(max_fun_eval/population_size);

% ### Optical Properties ###
opt_st.N_a_opt = sum(opt_st.face_alpha_opt); % number of aplha to optimize
opt_st.N_ei_opt = sum(opt_st.face_eps_int_opt); % number of eps_int to optimize
opt_st.N_ee_opt = sum(opt_st.face_eps_ext_opt); % number of eps_ext to optimize

% ### Conductance Links ###
% N_link X coord, N_link Y coord, N_link Z coord of the "from" list lower bounds (lb) and upper bounds (ub)
node_from = sat.node.globe(ismember([sat.node.globe.ID_item], find(opt_st.item_from_link))); % select nodes of the item_from_link
opt_st.id_node_from_list = [node_from.ID];
coord_from = [node_from.node];
x_from = coord_from(1:3:end);
y_from = coord_from(2:3:end);
z_from = coord_from(3:3:end);
coord_link_from_lb = [min(x_from), min(y_from), min(z_from)]; % x, y, z
coord_link_from_ub = [max(x_from), max(y_from), max(z_from)]; % x, y, z
% N_link X coord, N_link Y coord, N_link Z coord of the "to" list lower bounds (lb) and upper bounds (ub)
node_to = sat.node.globe(ismember([sat.node.globe.ID_item], find(opt_st.item_to_link))); % select nodes of the item_to_link
opt_st.id_node_to_list = [node_to.ID];
coord_to = [node_to.node];
x_to = coord_to(1:3:end);
y_to = coord_to(2:3:end);
z_to = coord_to(3:3:end);
coord_link_to_lb = [min(x_to), min(y_to), min(z_to)]; % x, y, z
coord_link_to_ub = [max(x_to), max(y_to), max(z_to)]; % x, y, z

% ### Heaters ###
opt_st.N_heater_hot = sum(opt_st.item_heater_hot == 1);
opt_st.N_heater_cold = sum(opt_st.item_heater_cold == 1);

% ### id of the design vector x ###

% N_a_opt + N_ei_opt + N_ee_opt + N_link*3 + N_link*3 + N_link*1 + N_heater_hot + N_heater_cold
% alpha     eps_int    eps_ext    coord_from coord_to   A_link     heater_hot_P   heater_cold_P 
opt_st.var_size = [opt_st.N_a_opt, opt_st.N_ei_opt, opt_st.N_ee_opt, opt_st.N_link*3, opt_st.N_link*3, opt_st.N_link*1, opt_st.N_heater_hot, opt_st.N_heater_cold];

opt_st.x_id_start = [];
opt_st.x_id_end = [];
it = 1;
for i = 1:1:length(opt_st.var_size)
opt_st.x_id_start = [opt_st.x_id_start, it];
it = it + opt_st.var_size(i);
opt_st.x_id_end = [opt_st.x_id_end, it - 1];
end

% ### lb and ub concatenation ###

% N_a_opt + N_ei_opt + N_ee_opt + N_link*3 + N_link*3 + N_link*1 + N_heater_hot + N_heater_cold
% alpha     eps_int    eps_ext    coord_from coord_to   A_link     heater_hot_P   heater_cold_P 

lb = [alpha_lb(opt_st.face_alpha_opt == 1),...
    eps_int_lb(opt_st.face_eps_int_opt == 1),...
    eps_ext_lb(opt_st.face_eps_ext_opt == 1),...
    repmat(coord_link_from_lb, 1, opt_st.N_link),...
    repmat(coord_link_to_lb, 1, opt_st.N_link),...
    repmat(A_link_lb, 1, opt_st.N_link),...
    heater_hot_lb,... (opt_st.item_heater_hot == 1) (opt_st.item_heater_cold == 1)
    heater_cold_lb];

ub = [alpha_ub(opt_st.face_alpha_opt == 1),...
    eps_int_ub(opt_st.face_eps_int_opt == 1),...
    eps_ext_ub(opt_st.face_eps_ext_opt == 1),...
    repmat(coord_link_from_ub, 1, opt_st.N_link),...
    repmat(coord_link_to_ub, 1, opt_st.N_link),...
    repmat(A_link_ub, 1, opt_st.N_link),...
    heater_hot_ub,... (opt_st.item_heater_hot == 1) (opt_st.item_heater_cold == 1)
    heater_cold_ub];

if sum(lb > ub) > 0
    error("Lower Bounds (lb) must be lower than Upper Bounds (ub)")
end

opt_st.N_var = opt_st.N_a_opt + opt_st.N_ei_opt + opt_st.N_ee_opt + opt_st.N_link*3 + opt_st.N_link*3 + opt_st.N_link*1 + opt_st.N_heater_hot + opt_st.N_heater_cold;
if opt_st.N_var == 0
    msgbox("Error: there are no variables for optimization.")
    error('Error: there are no variables for optimization.')
end

opt_st.obj_name = [];
opt_st.objective_active = zeros(1,5);
if val_cases == 2
    opt_st.N_objectives = 2;
    opt_st.obj_name = [opt_st.obj_name; "Temperature - Hot case"];
    opt_st.obj_name = [opt_st.obj_name; "Temperature - Cold case"];
    opt_st.objective_active(1:2) = 1;
else
    opt_st.N_objectives = 1;
    opt_st.obj_name = [opt_st.obj_name; "Temperature"];
    opt_st.objective_active(1) = 1;
end
if opt_st.N_heater_hot > 0
    opt_st.N_objectives = opt_st.N_objectives + 1;
    if val_cases == 2
        opt_st.obj_name = [opt_st.obj_name; "Heater power -  Hot case"];
    else       
        opt_st.obj_name = [opt_st.obj_name; "Heater power"]; 
    end
    opt_st.objective_active(3) = 1;
end
if opt_st.N_heater_cold > 0
    opt_st.N_objectives = opt_st.N_objectives + 1;
    opt_st.obj_name = [opt_st.obj_name; "Heater power -  Cold case"];
    opt_st.objective_active(4) = 1;
end
if opt_st.N_link > 0
    opt_st.N_objectives = opt_st.N_objectives + 1;
    opt_st.obj_name = [opt_st.obj_name; "Conductance link volume"];
    opt_st.objective_active(5) = 1;
end

%% Setup Var

% ### Operative Temperature Ranges and Weights ###
T_op_min_item = [];
T_op_max_item = [];
if isfield(sat.prop,"ext")
    T_op_min_item = [T_op_min_item, [sat.prop.ext.T_op_min]];
    T_op_max_item = [T_op_max_item, [sat.prop.ext.T_op_max]];
end
if isfield(sat.prop,"sp")
    T_op_min_item = [T_op_min_item, [sat.prop.sp.T_op_min]];
    T_op_max_item = [T_op_max_item, [sat.prop.sp.T_op_max]];
end
if isfield(sat.prop,"board")
    T_op_min_item = [T_op_min_item, [sat.prop.board.T_op_min]];
    T_op_max_item = [T_op_max_item, [sat.prop.board.T_op_max]];
end
if isfield(sat.prop,"parall")
    T_op_min_item = [T_op_min_item, [sat.prop.parall.T_op_min]];
    T_op_max_item = [T_op_max_item, [sat.prop.parall.T_op_max]];
end
if isfield(sat.prop,"cyl")
    T_op_min_item = [T_op_min_item, [sat.prop.cyl.T_op_min]];
    T_op_max_item = [T_op_max_item, [sat.prop.cyl.T_op_max]];
end
% convert to column vector
if size(T_op_max_item,2) > size(T_op_max_item,1)
    T_op_max_item = T_op_max_item';
end
if size(T_op_min_item,2) > size(T_op_min_item,1)
    T_op_min_item = T_op_min_item';
end
opt_st.T_op_item_avg = (T_op_min_item + T_op_max_item)/2; % average
opt_st.T_op_item = [T_op_min_item, T_op_max_item]; % min - max
range=(T_op_max_item-T_op_min_item);
opt_st.f=fit([min(range) max(range)]', [1 0.1]', 'poly1');
k_part=opt_st.f(range);
opt_st.k_weight=k_part/sum(k_part);

%% Delunay Search Preparation

% Conductance link, are thermal link from an internal node of each of the
% item in the "from" list to a node of one of the item in the "to" list.
% Conductance link may represent thermal straps or other thermal links.
% The optimization algorithms choose a set of (x,y,z) coord for the
% starting and finishing point of each strap; coordinates are converted to
% the nearest node id using Delunay triangulation search "dsearchn".

if opt_st.N_link > 0
    disp('Precomputation of the Delunay triangulation...');
    fictitious_node = [1e6, 1e6, 1e6]; % fictitious node to prevent degeneracy (coplanar points)
    % Triangulation of the "from" nodes
    opt_st.dsearchn_points_from = [[x_from', y_from', z_from']; fictitious_node];
    opt_st.Triang_from = delaunayn(opt_st.dsearchn_points_from);
    % Triangulation of the "to" nodes
    opt_st.dsearchn_points_to = [[x_to', y_to', z_to']; fictitious_node];
    opt_st.Triang_to = delaunayn(opt_st.dsearchn_points_to);
end

end