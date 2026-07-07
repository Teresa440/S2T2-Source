function [Source_term,sat]=Q_source3_compute(Source_term_partial,sat,r_dc,cases)

% Warning, a strong hypothesis is made to speed up computation:
% Hypothesis: external nodes have id ranging from 1 to N_nodes_ext

N_nodes = sat.node.total_node;
N_ext_node = sum([sat.node.ext.n_node]); % structure
if sat.geom.Nsp > 0
    N_ext_node = N_ext_node + sum([sat.node.sp.n_node]); % solar panels
end

t0 = length(r_dc(:,1)); % number of temporal instants
% return values
epsAf_earth = zeros(N_nodes, t0);
epsAf_space = zeros(N_nodes, t0);
Qs = zeros(t0, N_nodes);
Qa = zeros(t0, N_nodes);
Qir = zeros(t0, N_nodes);
% optical properties, exluding internal nodes
alpha_node_ext = zeros(t0, N_ext_node, Source_term_partial.max_faces_per_node);
eps_ext_node_ext = zeros(t0, N_ext_node, Source_term_partial.max_faces_per_node);

Source_term = struct;

for j=1:N_ext_node
    for k=1:length(sat.node.globe(j).face)
        face_id = sat.node.globe(j).face(k);
        % prop_opt 1,2,3 -> alpha, (eps_int), eps:
        alpha_node_ext(:,j,k) = sat.geom.surfaces(face_id).prop_opt(1);
        if isequal(sat.geom.surfaces(face_id).item, 'ex')
            eps_ext_node_ext(:,j,k) = sat.geom.surfaces(face_id).prop_opt(3); % structure has eps_ext in 3rd position
        elseif isequal(sat.geom.surfaces(face_id).item, 'sol')
            eps_ext_node_ext(:,j,k) = sat.geom.surfaces(face_id).prop_opt(2); % solar panels have eps_ext in 3rd position
        end
    end
end

%% Qir
M4 = Source_term_partial.Qir_partial_ext.*eps_ext_node_ext;                % size: t0*N_nodes_ext*max_faces_per_node
Qir_node_ext = sum(M4,3);                                                  % size: t0*N_nodes_ext

%% epsAf_earth
M4 = Source_term_partial.epsAf_earth_partial_ext.*eps_ext_node_ext;
epsAf_earth_node_ext = sum(M4,3);
epsAf_earth_node_ext = epsAf_earth_node_ext'; % transposition to match legacy data format

%% epsAf_space
M4 = Source_term_partial.epsAf_space_partial_ext.*eps_ext_node_ext;
epsAf_space_node_ext = sum(M4,3);
epsAf_space_node_ext = epsAf_space_node_ext'; % transposition to match legacy data format

%% Qa
M4 = Source_term_partial.Qa_partial_ext.*alpha_node_ext;
Qa_node_ext = sum(M4,3);

%% Qs
M4 = Source_term_partial.Qs_partial_ext.*alpha_node_ext;
Qs_node_ext = sum(M4,3);

%%
Qir(:,1:N_ext_node) = Qir_node_ext;
Qa(:,1:N_ext_node)= Qa_node_ext;
Qs(:,1:N_ext_node) = Qs_node_ext;
epsAf_earth(1:N_ext_node,:) = epsAf_earth_node_ext;
epsAf_space(1:N_ext_node,:) = epsAf_space_node_ext;
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