function [face] = flat_cap_face(rot,Center,is_top)

% Tappo piatto: due facce, normali costanti (nessuna suddivisione per
% anello/settore, non serve -- un disco piatto ha una sola normale per
% lato). face(1) = disco esterno (layer+core, verso 'ex'); face(2) = lato
% interno del core, verso la cavita' del cilindro (opposta a face(1)).
% is_top segue la stessa convenzione di build_cyl_cap.m: il tappo si
% estende in +z da z_ref per is_top=true (normale uscente = +z), in -z
% per is_top=false (normale uscente = -z), prima della rotazione.

face = struct('ID',{1,2},'norm',{[],[]},'mesh',{Center,Center},'gridX',{[],[]},'gridY',{[],[]},'gridZ',{[],[]});

if is_top
    face(1).norm = [0 0 1]*rot;
else
    face(1).norm = [0 0 -1]*rot;
end
face(2).norm = -face(1).norm;

end
