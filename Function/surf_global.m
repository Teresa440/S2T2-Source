function [sat,surf_for_MCRT_int,surf_for_MCRT_ext] = surf_global(sat)

ind=cellfun('isempty',{sat.node.globe.face});

facesss=sat.node.globe;
facesss(ind)=[];

n=length(sat.geom.globe);
surfaces=struct('ID',cell(1,n),...
            'id_node',cell(1,n),...
            'elem',cell(1,n),...
            'vert',cell(1,n),...
            'center',cell(1,n),...
            'match',cell(1,n),...
            'item',cell(1,n),...
            'number',cell(1,n),...
            'ex_in',cell(1,n),...
            'prop_opt',cell(1,n));

for i=1:1:n
    surfaces(i).ID=i;
end
for i=1:1:length(facesss)
    for j=1:1:length(facesss(i).face)
        id_f=facesss(i).face(j);
        surfaces(id_f).id_node=[surfaces(id_f).id_node,facesss(i).ID];
    end
end
for i=1:1:length(surfaces)
    surfaces(i).elem=sat.node.globe(surfaces(i).id_node);
end

for i=1:1:length(surfaces)
    for j=1:1:length(surfaces(i).elem)
        fx=length(surfaces(i).elem(j).face);
        if fx>=2
            for a=1:1:fx
            switch a
                case 1
                    if surfaces(i).elem(j).face(a)==i
                        surfaces(i).elem(j).vertf=surfaces(i).elem(j).vertf(1:4,:);                             
                        break
                    end
                case 2
                    if surfaces(i).elem(j).face(a)==i
                        surfaces(i).elem(j).vertf=surfaces(i).elem(j).vertf(5:8,:);                             
                        break
                    end                        
                case 3
                    if surfaces(i).elem(j).face(a)==i
                        surfaces(i).elem(j).vertf=surfaces(i).elem(j).vertf(9:12,:);                             
                        break
                    end                        
             end
            end
        end
    end
end


for i=1:1:length(surfaces)    
    surfaces(i).item=surfaces(i).elem(1).item;
    surfaces(i).number=surfaces(i).elem(1).number;
    surfaces(i).ex_in=surfaces(i).elem(1).ex_in;
    surfaces(i).norm=sat.geom.globe(i).norm;
    surfaces(i).prop_opt=sat.geom.globe(i).prop_opt;
end


for i=1:1:length(surfaces)
    cont=0;
    if strcmp(surfaces(i).item,'cyl')==0 && strcmp(surfaces(i).item,'sphcap')==0
        for j=1:1:length(surfaces(i).elem)
            if strcmp(surfaces(i).elem(j).type,'v')==1
                surfaces(i).vert=[surfaces(i).vert;surfaces(i).elem(j).node];
            end
        end
    elseif strcmp(surfaces(i).item,'cyl')==1
        % Gruppo 'cyl': raccolgo il vertf INTERO (non solo 2 punti) di
        % ogni elemento 's' (angolo), poi separo per raggio (distanza dal
        % baricentro) in anello esterno/interno. Un tappo piatto anulare
        % (R_int..R_out) e' un buco vero, non rappresentabile con un solo
        % contorno -- prima si prendevano solo 2 punti per elemento (il
        % "lato piu' esterno" della cella d'angolo), che per l'angolo
        % INTERNO (aggiunto quando e' stato esposto il foro della parete)
        % da' per errore il confine tra il 1' e 2' anello radiale (es.
        % r=70) invece del vero R_int (r=60): il poligono aggregato
        % risultava un contorno auto-intersecante che non escludeva il
        % foro (bug preesistente, gia' presente col solo angolo esterno,
        % un dodecagono a R_out che non ha mai escluso il centro) E
        % arrivava a escludere per errore parte della corona vera.
        % Verificato: fino all'83% dei raggi dalla calotta verso questo
        % gruppo si "impegnavano" su di esso (poligono aggregato troppo
        % permissivo verso il foro) per poi fallire il test per-elemento
        % e andare persi, invece di essere lasciati liberi per il vero
        % bersaglio.
        cand_full=[]; cand_simple=[]; n_s=0;
        for j=1:1:length(surfaces(i).elem)
            if strcmp(surfaces(i).elem(j).type,'s')==1
                cand_full=[cand_full;surfaces(i).elem(j).vertf]; %#ok
                cand_simple=[cand_simple;surfaces(i).elem(j).vertf(1:2,:)]; %#ok
                n_s=n_s+1;
            end
        end
        if isempty(cand_full)
            surfaces(i).vert=[];
        else
            % Un tappo piatto anulare aggrega SEMPRE gli elementi
            % d'angolo di TUTTI i settori (uno per i=1..Nt, quindi n_s
            % grande) -- una fascia LATERALE (un solo settore) ne ha al
            % piu' 2 (i due angoli sopra/sotto), 3 nel caso limite Nr==1.
            % Questo distingue in modo affidabile i due casi senza dover
            % conoscere Nt esplicitamente (surf_global.m non lo riceve).
            % Nota: provato prima un criterio puramente geometrico
            % (variazione di raggio perpendicolare alla normale) ma non
            % distingue i due casi -- per una fascia laterale il raggio e'
            % ANCHE li' costante "lungo la normale", stesso pattern
            % apparente di un tappo piatto; la vera differenza e'
            % strutturale (quanti settori contribuiscono), non geometrica.
            if n_s <= 4
                surfaces(i).vert = cand_simple;
            else
                axis_pt = mean(cand_full,1);
                delta = cand_full-axis_pt;
                along_normal = delta*surfaces(i).norm';
                perp = delta - along_normal*surfaces(i).norm;
                radii = vecnorm(perp,2,2);
                r_min = min(radii); r_max = max(radii);
            if (r_max-r_min) < 0.05*max(r_max,1)
                % variazione di raggio trascurabile: comportamento
                % originale invariato, un solo contorno.
                surfaces(i).vert = cand_simple;
            else
                % tappo anulare: due contorni separati, ciascuno solo i
                % punti del proprio bordo VERO (non i bordi interni delle
                % singole celle d'angolo), ordinati per angolo.
                tol = (r_max-r_min)*0.15;
                outer_pts = cand_full(radii >= r_max-tol,:);
                inner_pts = cand_full(radii <= r_min+tol,:);
                if size(outer_pts,1)>=3 && size(inner_pts,1)>=3
                    [~,outer_sorted] = center_sort_polygon(outer_pts);
                    [~,inner_sorted] = center_sort_polygon(inner_pts);
                    % NaN separa i due contorni: inpolygon() di MATLAB
                    % supporta nativamente questo formato per un poligono
                    % con un buco (dentro l'esterno, fuori dall'interno).
                    surfaces(i).vert = [outer_sorted; NaN(1,3); inner_sorted];
                else
                    surfaces(i).vert = cand_simple; % fallback prudente
                end
            end
            end
        end
    else
        for j=1:1:length(surfaces(i).elem)
            if strcmp(surfaces(i).elem(j).type,'s')==1 && cont==0
                surfaces(i).vert=[surfaces(i).vert;surfaces(i).elem(j).vertf(1:2,:)];
%                 surfaces(i).vert=[surfaces(i).vert;surfaces(i).elem(j).vertf(2:3,:)];


            elseif strcmp(surfaces(i).elem(j).type,'s')==1 && cont>=2
                surfaces(i).vert=[surfaces(i).vert;surfaces(i).elem(j).vertf(3:4,:)];
% %                 surfaces(i).vert=[surfaces(i).vert;surfaces(i).elem(j).vertf([1,4],:)];
            elseif strcmp(surfaces(i).elem(j).type,'cq')==1 && strcmp(surfaces(i).item,'sphcap')==1
                % Con la suddivisione fine per anello/blocco (calotta
                % sferica), la maggior parte dei gruppi non contiene
                % nessun elemento 's' (angolo) -- solo 'cq'. Senza questo
                % ramo il poligono aggregato resterebbe vuoto. Non tocca
                % 'cyl' (la parete ha sempre almeno un elemento 's' per
                % gruppo, comportamento invariato).
                %
                % Presi TUTTI e 4 i vertici, non solo 2: con
                % phi_block_size=1 (il valore attualmente usato in
                % GMM4.m/build_sphcap.m) ogni gruppo 'cq' ha esattamente UN
                % elemento, quindi il poligono del gruppo deve coincidere
                % col vero quadrilatero dell'elemento -- prenderne solo 2
                % vertici lo rendeva degenere (retta, area=0), facendo
                % fallire quasi sempre il pre-filtro inpolygon in
                % MC_ray_tracing3.m (fattore di vista della calotta
                % crollato fino al 50%, verificato). Se in futuro
                % phi_block_size>1 (piu' elementi 'cq' per gruppo), questo
                % prende 4 vertici per OGNI elemento invece di 2 --
                % ridondante ma non testato: ricontrollare che
                % l'ordinamento angolare (center_sort_polygon, sotto)
                % produca comunque il perimetro giusto in quel caso.
                surfaces(i).vert=[surfaces(i).vert;surfaces(i).elem(j).vertf];
            elseif (strcmp(surfaces(i).elem(j).type,'cb')==1 || strcmp(surfaces(i).elem(j).type,'ct')==1) ...
                    && strcmp(surfaces(i).item,'sphcap')==1
                % Gruppo dell'apice: un solo elemento ('cb'/'ct'), il cui
                % vertf e' gia' l'intero poligono di bordo (gli Nt punti
                % attorno al polo, vedi node_sphcap_creator.m), non un
                % lato condiviso con altri elementi -- va preso per
                % intero, non solo 2 punti come per 'cq'/'s'.
                surfaces(i).vert=[surfaces(i).vert;surfaces(i).elem(j).vertf];
            end
        end
          cont=cont+1;
    end
end



for i=1:1:length(surfaces)
     xyz=surfaces(i).vert;
     nan_row = find(isnan(xyz(:,1)),1);
     if ~isempty(nan_row)
         % Tappo anulare (guscio esterno/interno separati da NaN, vedi
         % sopra): centro/ordinamento dal solo contorno ESTERNO (quello
         % che definisce l'estensione reale della superficie); area =
         % esterno meno il foro. Il formato NaN-separato di .vert viene
         % preservato -- serve cosi' a MC_ray_tracing3.m per il test
         % inpolygon con buco.
         outer_part = xyz(1:nan_row-1,:);
         inner_part = xyz(nan_row+1:end,:);
         [xyzc,outer_sorted] = center_sort_polygon(outer_part);
         [~,inner_sorted] = center_sort_polygon(inner_part);
         surfaces(i).center = xyzc;
         surfaces(i).vert = [outer_sorted; NaN(1,3); inner_sorted];
         area_outer = area_polygon2(xyzc,outer_sorted(:,1),outer_sorted(:,2),outer_sorted(:,3));
         inner_c = mean(inner_sorted,1);
         area_inner = area_polygon2(inner_c,inner_sorted(:,1),inner_sorted(:,2),inner_sorted(:,3));
         surfaces(i).area = area_outer - area_inner;
     else
%      normal=surfaces(i).norm;
%      xyzc=mean(points,1);
    [xyzc,xyz] = center_sort_polygon(xyz);
    surfaces(i).center=xyzc;
   surfaces(i).vert=xyz;
   surfaces(i).vert= unique(surfaces(i).vert,'stable','rows');
%      [points] = sort_vert(points,normal);

%      surfaces(i).vert=points;
   surfaces(i).area=area_polygon2(xyzc,surfaces(i).vert(:,1),...
       surfaces(i).vert(:,2),surfaces(i).vert(:,3));
%    surfaces(i).normal=normal_from_points(surfaces(i).vert);
     end
end



for i=1:1:length(surfaces)
    for j=1:1:length(surfaces(i).elem)  
        xyz=surfaces(i).elem(j).vertf;
       [xyzc,xyz] = center_sort_polygon(xyz);
       surfaces(i).elem(j).center=xyzc;
       surfaces(i).elem(j).vertf=xyz;
       surfaces(i).elem(j).vertf= unique(surfaces(i).elem(j).vertf,'stable','rows');
       surfaces(i).elem(j).area=area_polygon2(xyzc,surfaces(i).elem(j).vertf(:,1),...
       surfaces(i).elem(j).vertf(:,2),surfaces(i).elem(j).vertf(:,3));
%        surfaces(i).elem(j).normal=normal_from_points(surfaces(i).elem(j).vertf);
    end
end



sat.geom.surfaces=surfaces;
surf_for_MCRT=surfaces;

%%
for s1=1:1:length(surf_for_MCRT)
    [Rot] = rot_mat_emiss(surf_for_MCRT(s1).norm, surf_for_MCRT(s1).center, surf_for_MCRT(s1).vert(1,:));
    surf_for_MCRT(s1).Rot = Rot;
    for i=1:1:length(surf_for_MCRT(s1).elem)
        coord_vert_2D = (Rot'*surf_for_MCRT(s1).elem(i).vertf')'; % result nx3
        coord_center_2D = (Rot'*surf_for_MCRT(s1).elem(i).center')'; % result 1x3
        x_0 = coord_center_2D(1);        
        y_0 = coord_center_2D(2); 
        diag = vecnorm(coord_vert_2D-coord_center_2D,2,2);
        R = max(diag);
        n_node_i = 8000;        
        r1 = rand(n_node_i,1);        
        r2 = rand(n_node_i,1);        
        r_i = R*sqrt(r1);        
        phi_i = 2*pi*r2;        
        x = x_0 + r_i.*cos(phi_i);        
        y = y_0 + r_i.*sin(phi_i);
        z = ones(n_node_i,1)*coord_center_2D(3);        
        flag = inpolygon(x,y,coord_vert_2D(:,1),coord_vert_2D(:,2));
        coord_i = [x(flag),y(flag),z(flag)];
        coord_i = Rot*coord_i';        
        surf_for_MCRT(s1).elem(i).p_start = coord_i'; 
    end
end

%%


surf_for_MCRT_int=surf_for_MCRT;
surf_for_MCRT_ext=surf_for_MCRT;

% EXCLUDE SOLAR PANEL FROM internal ANALYSIS :
ind=cellfun(@(v)any(strcmp(v,'sol')==1),{surf_for_MCRT_int.item});
surf_for_MCRT_int(ind)=[];

for i=1:1:length(surf_for_MCRT_int)
    surf_for_MCRT_int(i).ID=i;
    if strcmp(surf_for_MCRT_int(i).item,'ex')==1
        surf_for_MCRT_int(i).norm=-surf_for_MCRT_int(i).norm;
    end
end

% Classificazione geometrica "affaccia sulla cavita' (concava) del proprio
% oggetto" vs "affaccia verso l'esterno (convesso)": confronta la normale
% di ogni superficie con la direzione dal baricentro del proprio oggetto
% (stesso item+number) al centro della superficie stessa. Se puntano in
% direzioni opposte (prodotto scalare negativo) la superficie e' concava
% rispetto al proprio oggetto (es. foro interno del cilindro, guscio
% interno della calotta, bordo) -- se concordi, e' convessa (parete
% esterna, guscio esterno, tappi). Serve per evitare di abilitare lo
% scambio radiativo tra una faccia interna e una esterna dello STESSO
% oggetto (fisicamente impossibile, c'e' materiale solido in mezzo) --
% vedi nota sotto sull'esclusione stesso-item-stesso-number.
all_items = {surf_for_MCRT_int.item};
all_numbers = [surf_for_MCRT_int.number];
is_cavity_facing = false(1,length(surf_for_MCRT_int));
for i=1:1:length(surf_for_MCRT_int)
    same_obj = strcmp(all_items,surf_for_MCRT_int(i).item) & (all_numbers==surf_for_MCRT_int(i).number);
    centers_same_obj = reshape([surf_for_MCRT_int(same_obj).center],3,[])';
    obj_centroid = mean(centers_same_obj,1);
    radial = surf_for_MCRT_int(i).center - obj_centroid;
    is_cavity_facing(i) = dot(surf_for_MCRT_int(i).norm,radial) < 0;
end

for i=1:1:length(surf_for_MCRT_int)
    item1=surf_for_MCRT_int(i).item;
    num1=surf_for_MCRT_int(i).number;
    f_id1=surf_for_MCRT_int(i).ID;
   if strcmp(item1,'ex')==1
        for j=1:1:length(surf_for_MCRT_int)
            if surf_for_MCRT_int(j).ID ~= surf_for_MCRT_int(i).ID
            surf_for_MCRT_int(i).match=[surf_for_MCRT_int(i).match,surf_for_MCRT_int(j).ID];
            end
        end
   else
        for j=1:1:length(surf_for_MCRT_int)
            item2=surf_for_MCRT_int(j).item;
            num2=surf_for_MCRT_int(j).number;
            f_id2=surf_for_MCRT_int(j).ID;
            if f_id2 ~= f_id1
                % Esclusione stesso-item-stesso-number: corretta per un
                % oggetto convesso (pannello, board, parallelepipedo -- non
                % vede mai se stesso, per costruzione geometrica). Per
                % 'cyl'/'sphcap' l'oggetto puo' avere sia facce concave
                % (foro interno, guscio interno, bordo) sia convesse
                % (parete esterna, guscio esterno, tappi): le prime SI
                % vedono tra loro anche se stesso oggetto/numero, le
                % seconde no (e una concava non vede MAI una convessa dello
                % stesso oggetto -- ci sarebbe materiale solido in mezzo).
                % Vedi is_cavity_facing sopra.
                %
                % Trovato e corretto (era la causa del problema segnalato
                % qui prima): la prima versione di questo fix escludeva in
                % base al solo item, sbloccando per errore ANCHE le coppie
                % concava-convessa dello stesso oggetto (es. foro interno
                % vs parete esterna) -- MC_ray_tracing3.m non fa un vero
                % test di occlusione, quindi quei raggi "rubavano" hit al
                % vero bersaglio attraverso il materiale solido, facendo
                % crollare la somma dei fattori di vista della cavita'
                % (osservato: 2-18% invece di ~100%).
                same_item_no_self_view = strcmp(item1,item2)==1 && ...
                    ~(is_cavity_facing(i) && is_cavity_facing(j));
                if same_item_no_self_view
                    if num1 ~= num2
                      surf_for_MCRT_int(i).match=[surf_for_MCRT_int(i).match,surf_for_MCRT_int(j).ID];
                    end
                else
                    surf_for_MCRT_int(i).match=[surf_for_MCRT_int(i).match,surf_for_MCRT_int(j).ID];
                end
            end
        end
    end
end

% EXCLUDE Boards Parall and Cyl FROM external ANALYSIS :
ind=cellfun(@(v)any(strcmp(v,'board')==1),{surf_for_MCRT_ext.item});
surf_for_MCRT_ext(ind)=[];
ind=cellfun(@(v)any(strcmp(v,'paral')==1),{surf_for_MCRT_ext.item});
surf_for_MCRT_ext(ind)=[];
ind=cellfun(@(v)any(strcmp(v,'cyl')==1),{surf_for_MCRT_ext.item});
surf_for_MCRT_ext(ind)=[];
ind=cellfun(@(v)any(strcmp(v,'sphcap')==1),{surf_for_MCRT_ext.item});
surf_for_MCRT_ext(ind)=[];

for i=1:1:length(surf_for_MCRT_ext)
    item1=surf_for_MCRT_ext(i).item;
    num1=surf_for_MCRT_ext(i).number;
    f_id1=surf_for_MCRT_ext(i).ID;  
    for j=1:1:length(surf_for_MCRT_ext)
        item2=surf_for_MCRT_ext(j).item;
        num2=surf_for_MCRT_ext(j).number;
        f_id2=surf_for_MCRT_ext(j).ID;
        if f_id2 ~= f_id1
            if strcmp(item1,item2)==1
                if num1 ~= num2
                  surf_for_MCRT_ext(i).match=[surf_for_MCRT_ext(i).match,surf_for_MCRT_ext(j).ID];
                end
            else
                surf_for_MCRT_ext(i).match=[surf_for_MCRT_ext(i).match,surf_for_MCRT_ext(j).ID];
            end
        end
    end
end













end