function [sat,eps_for_Ge] = opt_prop(sat)

nn = sat.node.total_node;
eps_for_Ge = zeros(nn,1);
for i = 1:1:nn
    sat.node.globe(i).n_face = length(sat.node.globe(i).face);
end
for i=1:1:length(sat.geom.surfaces)
    node_f = sat.geom.surfaces(i).id_node;
    % for the structure corresponds to eps int, for the other item to the only eps
    eps_for_Ge(node_f) = eps_for_Ge(node_f) + sat.geom.surfaces(i).prop_opt(2);
end
id_prop_eps = eps_for_Ge ~= 0; % every node that has at least 1 surface
n_face = [sat.node.globe.n_face]';
eps_for_Ge(id_prop_eps) = eps_for_Ge(id_prop_eps)./n_face(id_prop_eps);

for i=1:1:nn
    if strcmp(sat.node.globe(i).item,'ex')==1
        sat.node.globe(i).Af_tot=sum(sat.node.globe(i).Af)/2; % only internal-facing area
    elseif strcmp(sat.node.globe(i).item,'sol')==1 || strcmp(sat.node.globe(i).item,'board')==1
        sat.node.globe(i).Af_tot=sum(sat.node.globe(i).Af)*2; % top and bottom sides
    else 
        sat.node.globe(i).Af_tot=sum(sat.node.globe(i).Af);
    end

    % sat.node.globe(i).prop_opt(2) = eps_for_Ge(i);
end

end