function [sat]=Q_dissipation3(sat,diss_mat,val_cases)

% Warning, strong hypothesis are made to speed up computation:
% Hypothesis: external nodes have id ranging from 1 to N_nodes_ext
% Hypothesis: sat.node.coord_nod is sorted by face id
% Hypothesis: and face id ranges from 1 to 6 + sat.NB

N_item = 1 + sat.geom.Nsp + sat.geom.Nb + sat.geom.NP + sat.geom.Nc;
Q_diss = zeros(sat.node.total_node, 1); % W of dissipation per node
node_to_ID_item = [sat.node.globe.ID_item];
V_nodes = [sat.node.globe.V];

if val_cases == 1 % hot/default case is second column (fisrt col is ID item)
    diss_mat = diss_mat(:,2);
elseif val_cases == 2 % cold case is second column (fisrt col is ID item)
    diss_mat = diss_mat(:,3);
end

for i = 1:1:N_item
    if ~isequal(diss_mat(i),0)
        node_diss = find(node_to_ID_item == i); % node that dissipates
        V_node_item = V_nodes(node_diss); % volumes of nodes of item i that dissipates
        Q_diss(node_diss)=diss_mat(i).*V_node_item./sum(V_node_item); % equal distribution with nodes
    end
end

if val_cases == 1
    sat.node.Q_diss=Q_diss;
elseif val_cases == 2
    sat.node.Q_diss_cold=Q_diss;
else
    disp("Error: val cases must be 1 or 2")
end

end