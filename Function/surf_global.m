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
    if strcmp(surfaces(i).item,'cyl')==0
        for j=1:1:length(surfaces(i).elem)
            if strcmp(surfaces(i).elem(j).type,'v')==1
                surfaces(i).vert=[surfaces(i).vert;surfaces(i).elem(j).node];
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
            end
        end 
          cont=cont+1;
    end   
end



for i=1:1:length(surfaces)
     xyz=surfaces(i).vert;
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
                if strcmp(item1,item2)==1
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