function [Vf,VV] = MC_ray_tracing2(surf_for_MCRT,sat,ButtonHandle)

surface=surf_for_MCRT;
tot_nod=length(sat.node.globe(1,:));
ns=length(surface);

Vf=zeros(tot_nod,tot_nod);
rr=zeros(tot_nod,1);
VV=struct('pend',cell(1,tot_nod));

rays=5000;

% CHANGE: added this loop
for s1=1:1:length(surface)
    [Rot] = rot_mat_emiss(surface(s1).norm, surface(s1).center, surface(s1).vert(1,:));
    surface(s1).Rot = Rot;
end

for s1=1:1:length(surface)

    for i=1:1:length(surface(s1).elem)
        % CHANGE: until "for r=1:1:rays"
%         norm1=surface(s1).norm;
%         p_start=surface(s1).elem(i).center;
%         ei=surface(s1).elem(i).ID;
%         p=surface(s1).elem(i).vertf(1,:);
%         [Rot] = rot_mat_emiss(norm1,p_start,p);
        p_start=surface(s1).elem(i).center;
        ei=surface(s1).elem(i).ID;
        Rot = surface(s1).Rot;
        for r=1:1:rays
            
            if strcmp(surface(s1).item,'board')==1
                theta=pi*rand(1);
                phi=2*pi*rand(1);
                omegaG=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];                
                omega=Rot*omegaG; % CHANGE
            else
                [omega] = random_dir(Rot); % omega=R*omegaG;
            end
            d_old=1000000;
            rr(ei)=rr(ei)+1;
            s2_t=0;
            for m=1:1:length(surface(s1).match)
                
                s2=surface(s1).match(m);
                vert_s=surface(s2).vert;
                n_plane=surface(s2).norm;
                p_plane=surface(s2).center;
                % area_s=surface(s2).area; % CHANGE 
                
                [p_end,d] = line_plane_inters(p_start,omega,n_plane,p_plane);
                % CHANGE: until "if"
                Rot_s2 = surface(s2).Rot;
                p_end_rot = (Rot_s2')*(p_end'); % 3x3 * 3x1 = 3x1
                vert_s_rot = (Rot_s2')*(vert_s'); % 3x3 * 3x4 = 3x4
                flag1 = inpolygon(p_end_rot(1),p_end_rot(2),vert_s_rot(1,:),vert_s_rot(2,:));
                %[flag1] = is_inside2(area_s,vert_s(:,1),vert_s(:,2),vert_s(:,3),p_end);
                if flag1==1 && d>0 && abs(d)<abs(d_old)
                    d_old=d;
                    s2_t=s2;
                    p_rec=p_end;
                end
            end
            
            if s2_t ~= 0
                % ref
%                 flag2 = 0;
                VV(ei).pend=[VV(ei).pend;p_rec];
                for j=1:1:length(surface(s2_t).elem)
                    vert_e=surface(s2_t).elem(j).vertf;
                    ej=surface(s2_t).elem(j).ID;
                    % area_e=surface(s2_t).elem(j).area; % CHANGE
                    % CHANGE: until "if"
                    Rot_s2_t = surface(s2_t).Rot;
                    p_rec_rot = (Rot_s2_t')*(p_rec'); % 3x3 * 3x1
                    vert_e_rot = (Rot_s2_t')*(vert_e'); % 3x3 * 3x4
                    flag2 = inpolygon(p_rec_rot(1),p_rec_rot(2),vert_e_rot(1,:),vert_e_rot(2,:));
                    %[flag2] = is_inside2(area_e,vert_e(:,1),vert_e(:,2),vert_e(:,3),p_rec);
                    if flag2==1
                        Vf(ei,ej)=Vf(ei,ej)+1;
                        break
                        
                    end
                end
%                 if flag2 == 0
%                     disp("NO elem of s2_t found")
%                 end
%             else
%                 disp("NO Intersection")
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