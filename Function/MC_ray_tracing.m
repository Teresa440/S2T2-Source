function [Vf,VV] = MC_ray_tracing(surf_for_MCRT,sat,ButtonHandle)

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
        p_start=surface(s1).elem(i).center;
        ei=surface(s1).elem(i).ID;
        p=surface(s1).elem(i).vertf(1,:);
        [Rot] = rot_mat_emiss(norm1,p_start,p);
        for r=1:1:rays
            %
            if strcmp(surface(s1).item,'board')==1 || strcmp(surface(s1).item,'sol')==1
                theta=pi*rand(1);
                phi=2*pi*rand(1);
                omega=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
            else
                [omega] = random_dir(Rot);
            end
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
                % ref
                VV(ei).pend=[VV(ei).pend;p_rec];
                for j=1:1:length(surface(s2_t).elem)
                    vert_e=surface(s2_t).elem(j).vertf;
                    ej=surface(s2_t).elem(j).ID;
                    area_e=surface(s2_t).elem(j).area;
                    [flag2] = is_inside2(area_e,vert_e(:,1),vert_e(:,2),vert_e(:,3),p_rec);
                    if flag2==1
                        Vf(ei,ej)=Vf(ei,ej)+1;
                        break
                        
                    end
                end
            end
%             else
%                 disp("NO Intersection")
%             end
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
        drawnow()
    end
    disp(perc*100)
end

for i=1:1:tot_nod
    if rr(i)>0
    Vf(i,:)=Vf(i,:)/rr(i);
    end
end

end