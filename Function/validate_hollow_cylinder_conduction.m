% Validazione: conduzione radiale in un cilindro cavo omogeneo, regime stazionario.
%
% Confronta la resistenza termica radiale "vista" dalla rete di conduttanze
% del modello (Circle_Mesh -> Mesh2D_to_Mesh3D -> node_cyl_creator3 -> TMM2)
% con la soluzione analitica esatta di un guscio cilindrico:
%
%   R_analitica = ln(r2/r1) / (2*pi*k*L)
%
% Nota: la rete di conduttanze collega solo i BARICENTRI degli anelli, non le
% superfici fisiche vere e proprie (r1=R_int e r2=R_out). Mancano quindi due
% "mezze resistenze" di bordo (dal baricentro del primo/ultimo anello alla
% superficie fisica), che nel modello sono accoppiate solo via radiazione
% (campo Af), non via conduzione. Per un confronto bordo-fisico-a-bordo-fisico
% coerente con la formula a mano, queste due mezze celle vengono aggiunte qui
% esplicitamente, con la stessa formula logaritmica usata in TMM2.m.

clear
addpath(fileparts(mfilename('fullpath')));

%% Parametri del caso di validazione (modificabili)
R_int = 50;    % raggio interno del foro [mm]
R_out = 100;   % raggio esterno [mm]
Lc    = 200;   % lunghezza del cilindro [mm]
k     = 200;   % conducibilita' termica [W/(m*K)]
Nt = 36;       % settori angolari
Nr = 20;       % anelli radiali
Nz = 2;        % livelli assiali (Nz=2 -> un solo strato, problema puramente radiale)

%% Costruzione della mesh e della rete di conduttanze
zz=linspace(-Lc/2,Lc/2,Nz);
[Nodes, Triangles, Quads]=Circle_Mesh(R_out,Nr,Nt,R_int);
[Nodes3D,~,Bricks] = Mesh2D_to_Mesh3D(Nodes,Triangles,Quads,zz);
Central=[];
total_nodes=size(Bricks,1);
[elem,Con] = node_cyl_creator3(Nodes3D,Central,Bricks,R_out,Lc,Nt,Nr,Nz,total_nodes,R_int);

for kk=1:total_nodes
    elem(kk).item='cyl';
    elem(kk).number=1;
    elem(kk).prop_mech=[2700e-9,900,k]; % [rho, cp, conducibilita']
    elem(kk).Af_tot=sum(elem(kk).Af);
end
sat.node.total_node=total_nodes;
sat.node.globe=elem;
sat.geom.cyl(1).Nt=Nt; sat.geom.cyl(1).L=Lc; sat.geom.cyl(1).Nz=Nz;
Vf_G=zeros(total_nodes); eps_int=zeros(1,total_nodes); sigma=5.67e-8;

[Gc,~,~] = TMM2(sat,sigma,Vf_G,Con,eps_int);

%% Riduzione della rete: percorso radiale puro, un solo settore (i=1), un solo strato assiale (h=1)
h=1; i=1;
kidx=@(i,j,h) Nt*Nr*(h-1) + (j-1)*Nt + i; % stessa formula di node_cyl_creator3 (caso cavo)

R_baricentro_a_baricentro = 0;
for j=1:Nr-1
    m1=kidx(i,j,h); m2=kidx(i,j+1,h);
    R_baricentro_a_baricentro = R_baricentro_a_baricentro + 1/Gc(m1,m2);
end

alfa=(360/Nt)/2; dz=Lc/(Nz-1); Cg=2*sind(alfa)*dz*1e-3; % [m]: A(r)=Cg*r per questa mesh

elem1  = elem(kidx(i,1,h));
elemNr = elem(kidx(i,Nr,h));
R_half_bordo_interno  = 0.5*log(elem1.Ac(2)/elem1.Af(5))  /(k*Cg); % Ac(5)=0 sul foro: usa Af
R_half_bordo_esterno  = 0.5*log(elemNr.Af(2)/elemNr.Ac(5))/(k*Cg); % Ac(2)=0 sul bordo est.: usa Af

R_settore_rete = R_half_bordo_interno + R_baricentro_a_baricentro + R_half_bordo_esterno;
R_cilindro_rete = R_settore_rete/Nt; % Nt settori identici in parallelo

%% Soluzione analitica esatta
R_cilindro_analitico = log(R_out/R_int)/(2*pi*k*(Lc*1e-3));

%% Report
fprintf('--- Validazione conduzione radiale, cilindro cavo ---\n');
fprintf('R_int=%g mm, R_out=%g mm, L=%g mm, k=%g W/mK, Nt=%d, Nr=%d\n\n', R_int,R_out,Lc,k,Nt,Nr);
fprintf('Resistenza cilindro intero (rete, baricentro-a-baricentro): %.6f K/W\n', R_baricentro_a_baricentro/Nt);
fprintf('Resistenza cilindro intero (rete, bordo-fisico-a-bordo):    %.6f K/W\n', R_cilindro_rete);
fprintf('Resistenza cilindro intero (analitica, ln(r2/r1)/(2 pi k L)): %.6f K/W\n', R_cilindro_analitico);
fprintf('Errore relativo (bordo-a-bordo vs analitica): %.4f%%\n', ...
    (R_cilindro_rete-R_cilindro_analitico)/R_cilindro_analitico*100);
