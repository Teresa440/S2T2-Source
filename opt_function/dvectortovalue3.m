function [alpha_str,eps_int_str,eps_ext_str,heat,id_node_from,id_node_to,A_link] = dvectortovalue3(x,opt_st)

n_sol = size(x,1); % number of solution in the pareto front
heat.case1 = zeros(opt_st.x_id_end(7) - opt_st.x_id_start(7) + 1, n_sol + 1); % warning, heat.case1/2 are saved transposed to match legacy data format
heat.case2 = zeros(opt_st.x_id_end(8) - opt_st.x_id_start(8) + 1, n_sol + 1); % warning, heat.case1/2 are saved transposed to match legacy data format
id_node_from = zeros(n_sol, opt_st.N_link);
id_node_to = zeros(n_sol, opt_st.N_link);

% deconstructing the design vector x
alpha_str = x(:, opt_st.x_id_start(1):opt_st.x_id_end(1));
eps_int_str = x(:, opt_st.x_id_start(2):opt_st.x_id_end(2));
eps_ext_str = x(:, opt_st.x_id_start(3):opt_st.x_id_end(3));
A_link = x(:, opt_st.x_id_start(6):opt_st.x_id_end(6));
%heater hot
heat.case1(:,1) = find(opt_st.item_heater_hot); % first column of heat is the id of the item
heat.case1(:,2:end) = x(:, opt_st.x_id_start(7):opt_st.x_id_end(7))'; % hot/default (note: transpose)
% heater cold
heat.case2(:,1) = find(opt_st.item_heater_cold); % first column of heat is the id of the item
heat.case2(:,2:end) = x(:, opt_st.x_id_start(8):opt_st.x_id_end(8))'; % cold (note: transpose)

if opt_st.N_link > 0
    for i = 1:1:n_sol
        x_one = x(i,:);
        coord_from = x_one(opt_st.x_id_start(4):opt_st.x_id_end(4)); % size = 1 x (N_link*3)
        coord_to = x_one(opt_st.x_id_start(5):opt_st.x_id_end(5));   % size = 1 x (N_link*3)
        % reshaping coords
        coord_from = reshape(coord_from, 3, opt_st.N_link)'; % size = N_link x 3; columns 1, 2, 3 are x, y, z of the points, respecively
        coord_to = reshape(coord_to, 3, opt_st.N_link)';     % size = N_link x 3; columns 1, 2, 3 are x, y, z of the points, respecively
        % id nodes of the "from" and "to" list
        id_node_from(i,:) = dsearchn(opt_st.dsearchn_points_from,opt_st.Triang_from,coord_from);
        id_node_to(i,:) = dsearchn(opt_st.dsearchn_points_to,opt_st.Triang_to,coord_to);
    end
else
    id_node_from = [];
    id_node_to = [];
end

end