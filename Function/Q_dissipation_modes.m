function [sat] = Q_dissipation_modes(sat,diss_mat,val_cases)

N_item = 1 + sat.geom.Nsp + sat.geom.Nb + sat.geom.NP + sat.geom.Nc;
node_to_ID_item = [sat.node.globe.ID_item];
V_nodes = [sat.node.globe.V];

if val_cases == 1 % hot/default case is second column (fisrt col is ID item)
   N_modes = size(diss_mat, 2) - 1;
   Q_diss_modes = zeros(sat.node.total_node,N_modes);
       for m = 1:N_modes
        col = m + 1;
        for i = 1:N_item
            if ~isequal(diss_mat(i, col), 0)
                node_diss = find(node_to_ID_item == i); %node that dissipates
                V_node_item = V_nodes(node_diss); %volumes of nodes of item i that dissipates
                Q_diss_modes(node_diss, m) = diss_mat(i, col) .* V_node_item ./ sum(V_node_item); % equal distribution with nodes
            end
        end
    end

    sat.node.Q_diss_modes = Q_diss_modes;


elseif val_cases == 2 % cold case is second column (fisrt col is ID item)
    N_modes = (size(diss_mat, 2) - 1) / 2;
    Q_diss_modes_hot  = zeros(sat.node.total_node, N_modes);
    Q_diss_modes_cold = zeros(sat.node.total_node, N_modes);

    for m = 1:N_modes
        col_hot  = 2*m;
        col_cold = 2*m + 1;
        for i = 1:N_item
            node_diss = find(node_to_ID_item == i);
            V_node_item = V_nodes(node_diss);

            if ~isequal(diss_mat(i, col_hot), 0)
                Q_diss_modes_hot(node_diss, m) = diss_mat(i, col_hot) .* V_node_item ./ sum(V_node_item);
            end
            if ~isequal(diss_mat(i, col_cold), 0)
                Q_diss_modes_cold(node_diss, m) = diss_mat(i, col_cold) .* V_node_item ./ sum(V_node_item);
            end
        end
    end

    sat.node.Q_diss_modes = Q_diss_modes_hot;
    sat.node.Q_diss_modes_cold = Q_diss_modes_cold;

else
    disp("Error: val cases must be 1 or 2")
end

end





