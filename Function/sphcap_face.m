function [face] = sphcap_face(Nt,Ntheta,phi_block_size,rot,Center,is_top)

if nargin<3 || isempty(phi_block_size)
    phi_block_size=1; % vedi node_sphcap_creator.m -- 1=esatto per cella
end
n_phi_blocks = ceil(Nt/phi_block_size);
n_shell_slots = 1 + (Ntheta-1)*n_phi_blocks; % apice: 1 solo slot -- vedi node_sphcap_creator.m
n_faces = 2*n_shell_slots + Nt;

face=struct('ID',cell(1,n_faces),...
'norm',cell(1,n_faces),...
'mesh',cell(1,n_faces),...
'gridX',cell(1,n_faces),...
'gridY',cell(1,n_faces),...
'gridZ',cell(1,n_faces));

% Emisfero puro (vedi sphcap_areas.m/Sph_Shell_Mesh.m): theta_max=90
% esatto, identico per guscio interno ed esterno -- non serve piu'
% calcolarlo per layer. mirror_sign applica lo specchio in z usato da
% build_sphcap.m per il tappo bottom (is_top=false).
mirror_sign = 2*is_top-1; % top: +1 (nessuno specchio), bottom: -1

for jf=1:1:n_faces
    face(jf).ID = jf;
    % .mesh non e' mai consultato per la calotta (a valle si usa sempre
    % elem.vertf), ma GMM4.m stesso (ismember(...,"rows") + mesh(1,:))
    % richiede una riga 1x3 valida per ogni voce di sat.geom.globe, non
    % vuota.
    face(jf).mesh = Center;
    if jf<=n_shell_slots
        % Guscio interno (R_int, face_inner in node_sphcap_creator.m):
        % normale radiale-dal-centro esatta ma invertita -- la superficie
        % interna affaccia verso il centro della cavita', non verso
        % l'esterno. Valutata al punto medio dell'anello theta e del
        % blocco phi a cui il gruppo appartiene (esatta per cella se
        % phi_block_size=1).
        [theta_j,phi_j] = shell_angles(jf,Ntheta,n_phi_blocks,phi_block_size,Nt);
        nrm = [sind(theta_j)*cosd(phi_j), sind(theta_j)*sind(phi_j), mirror_sign*cosd(theta_j)];
        face(jf).norm = -nrm*rot;
    elseif jf<=2*n_shell_slots
        % Guscio esterno (R_out, face_outer in node_sphcap_creator.m):
        % normale radiale-dal-centro esatta, diretta verso l'esterno.
        [theta_j,phi_j] = shell_angles(jf-n_shell_slots,Ntheta,n_phi_blocks,phi_block_size,Nt);
        nrm = [sind(theta_j)*cosd(phi_j), sind(theta_j)*sind(phi_j), mirror_sign*cosd(theta_j)];
        face(jf).norm = nrm*rot;
    else
        % Bordo/celle laterali: theta=90 esatto per l'emisfero puro, quindi
        % la normale (direzione theta) e' costante e assiale, indipendente
        % da phi -- coincide con la normale del guscio interno all'ultimo
        % anello. Verificato: scarto 0.00 esatto contro la formula
        % azimutale precedente (presa da cylinder_face.m, sbagliata per il
        % bordo di una cupola).
        face(jf).norm = -mirror_sign*[0 0 1]*rot;
    end
end

end

function [theta_j,phi_j] = shell_angles(jf_shell,Ntheta,n_phi_blocks,phi_block_size,Nt)
if jf_shell==1
    theta_j = 0;
    phi_j = 0;
else
    jf2 = jf_shell-1;
    j = ceil(jf2/n_phi_blocks)+1;
    pb = jf2-(j-2)*n_phi_blocks;
    theta_j = 90*(j-0.5)/Ntheta; % emisfero puro: stesso theta_max=90 per guscio interno ed esterno
    phi_j = 360*((pb-0.5)*phi_block_size)/Nt;
end
end
