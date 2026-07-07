function [Vf,VV,VV_ref] = MC_RT_Reflection(surf_for_MCRT,sat,ButtonHandle)

surface=surf_for_MCRT;
tot_nod=length(sat.node.globe(1,:));
ns=length(surface);

Vf=zeros(tot_nod,tot_nod);
rr=zeros(tot_nod,1);
VV=struct('pend',cell(1,tot_nod));

rays=5000;

for s1=1:1:length(surface)
% 
    for i=1:1:length(surface(s1).elem)
        %
        norm1=surface(s1).norm;
%         p_start=surface(s1).elem(i).center;
        ei=surface(s1).elem(i).ID;
        p=surface(s1).elem(i).vertf(1,:);
        [Rot] = rot_mat_emiss(norm1,surface(s1).elem(i).center,p);
        for r=1:1:rays            
            id_p_start = randi(length(surface(s1).elem(i).p_start));
            p_start=surface(s1).elem(i).p_start(id_p_start,:);
            % [flag4] = is_inside2(surface(s1).elem(i).area,...
            %     surface(s1).elem(i).vertf(:,1),surface(s1).elem(i).vertf(:,2),surface(s1).elem(i).vertf(:,3),p_start);
            if strcmp(surface(s1).item,'board')==1 || strcmp(surface(s1).item,'sol')==1
                theta=acos(2*rand(1)-1);
                phi=2*pi*rand(1);
                omega=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
            else
                [omega] = random_dir(Rot);
               
            end
%             Om=p_start+omega'*10;
%             plot3([p_start(:,1);Om(:,1)],[p_start(:,2);Om(:,2)],[p_start(:,3);Om(:,3)],'linewidth',5)
            d_old=1000000;
            rr(ei)=rr(ei)+1;

            s2_t=0;
            for m=1:1:length(surface(s1).match)
                %
                s2=surface(s1).match(m);
                vert_s=surface(s2).vert;
                n_plane=surface(s2).norm;
                p_plane=surface(s2).center;
                area_s=surface(s2).area;  
                %
                [p_end,d] = line_plane_inters(p_start,omega,n_plane,p_plane);
                [flag1] = is_inside2(area_s,vert_s(:,1),vert_s(:,2),vert_s(:,3),p_end);
                if flag1==1 && d>0 && abs(d)<abs(d_old)
                    d_old=d;
                    s2_t=s2;
                    p_rec=p_end;
                end
                
            end
            %
            if s2_t ~= 0
                VV(ei).pend=[VV(ei).pend;p_rec];
                %reflection
                Ra=rand;
                s_end=0;
                if surface(s2_t).prop_opt(2) < Ra
                    s_start=s2_t;
                    p_start_ref=p_rec;
                    [s_end,p_end_ref,n_end,ref,flag3] = reflection(surf_for_MCRT,s_start,p_start_ref,p_start);
                    p_rec=p_end_ref;
                    if flag3==1
                        r=r-1;
                        rr(ei)=rr(ei)-1;
                        continue
                    end
                    VV_ref(ei).Ray(rays).ref=ref;
                end
                if s_end > 0
                   s2_t=s_end;
                end
                contato=0;
                for j=1:1:length(surface(s2_t).elem)
                    vert_e=surface(s2_t).elem(j).vertf;
                    ej=surface(s2_t).elem(j).ID;
                    area_e=surface(s2_t).elem(j).area;
                    [flag2] = is_inside2(area_e,vert_e(:,1),vert_e(:,2),vert_e(:,3),p_rec);
                    
                    if flag2==1
                        Vf(ei,ej)=Vf(ei,ej)+1;
                        contato=1;
                         break
                        
                    end
                end

            end              
            
            
        end

    end
    perc=s1/ns;
    if exist('ButtonHandle','var')        
        ButtonHandle.Text = "Surface " + string(s1) + "/" + string(ns);
        currentProg = min(round((size(ButtonHandle.Icon,2)-2)*(perc)),size(ButtonHandle.Icon,2)-2);    
        RGB = ButtonHandle.Icon;    
        RGB(2:end-1, 2:currentProg+1, 1) = 6/255;
        RGB(2:end-1, 2:currentProg+1, 2) = 176/255;
        RGB(2:end-1, 2:currentProg+1, 3) = 37/255;
        ButtonHandle.Icon = RGB;
        drawnow
    end
    disp(perc*100)
end

for i=1:1:tot_nod
    if rr(i)>0
    Vf(i,:)=Vf(i,:)/rr(i);
    end
end

end