clear
addpath(fileparts(mfilename('fullpath')));

%% Parametri
R_int = 50; R_out = 100; Lc = 200; k = 200;
Nt = 36; Nr = 20; Nz = 2;

%% Costruzione mesh e rete di conduttanze (stessa pipeline del tool)
zz = linspace(-Lc/2, Lc/2, Nz);
[Nodes, Triangles, Quads] = Circle_Mesh(R_out, Nr, Nt, R_int);
[Nodes3D, ~, Bricks] = Mesh2D_to_Mesh3D(Nodes, Triangles, Quads, zz);
total_nodes = size(Bricks,1);
[elem, Con] = node_cyl_creator3(Nodes3D, [], Bricks, R_out, Lc, Nt, Nr, Nz, total_nodes, R_int);
total_nodes = numel(elem); % may include 2 extra sealing-cap nodes

for kk = 1:total_nodes
    elem(kk).item = 'cyl';
    elem(kk).number = 1;
    elem(kk).prop_mech = [2700, 900, k];
    elem(kk).Af_tot = sum(elem(kk).Af);
end
sat.node.total_node = total_nodes;
sat.node.globe = elem;
sat.geom.cyl(1).Nt = Nt; sat.geom.cyl(1).L = Lc; sat.geom.cyl(1).Nz = Nz;
Vf_G = zeros(total_nodes); eps_int = zeros(1,total_nodes); sigma = 5.67e-8;

[Gc,~,~] = TMM2(sat, sigma, Vf_G, Con, eps_int);

%% Estrai la catena radiale (un solo settore, i=1, unico strato h=1)
kidx = @(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i;
i = 1; h = 1;

%% Resistenze dei singoli collegamenti e resistenza cumulata
R_link = zeros(1,Nr-1);
for j = 1:Nr-1
    m1 = kidx(i,j,h); m2 = kidx(i,j+1,h);
    R_link(j) = 1/Gc(m1,m2);
end
R_cum = [0, cumsum(R_link)];   % R_cum(j) = resistenza dal nodo 1 al nodo j

%% Temperature del modello (partitore di resistenze)
T_hot = 100; T_cold = 0;
T_model = T_hot - (T_hot-T_cold) * R_cum / R_cum(end);

%% Raggi efficaci di ogni anello (calcolati indipendentemente, a mano)
a = (R_out - R_int)/Nr;
r_in  = R_int + (0:Nr-1)*a;
r_out = R_int + (1:Nr)*a;
r_eff = sqrt(r_in .* r_out);      % media geometrica: dove "vive" fisicamente ogni nodo

%% Temperature analitiche, valutate agli stessi raggi efficaci
T_analytic = T_hot - (T_hot-T_cold) * log(r_eff/r_eff(1)) / log(r_eff(end)/r_eff(1));

%% Confronto
err = abs(T_model - T_analytic);
fprintf('Errore massimo assoluto: %.6f °C\n', max(err));
fprintf('Errore massimo relativo: %.6f %%\n', max(err)/(T_hot-T_cold)*100);

figure
plot(r_eff, T_model, 'o-', r_eff, T_analytic, 'x--')
legend('Modello (rete)', 'Analitica (log)')
xlabel('raggio [mm]'); ylabel('T [°C]')
title('Profilo di temperatura radiale: modello vs analitica')