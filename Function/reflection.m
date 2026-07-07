function [s_end,p_end,n_end,ref,flag3] = reflection(surf_for_MCRT,s_start,p_start,p_dir)
surface=surf_for_MCRT;
s=s_start;
flag3=0;
if strcmp(surface(s).item,'board')==1
    retta=(p_dir-p_start)/norm(p_dir-p_start);
    if dot(retta,surface(s).norm) < 0
        surface(s).norm=-surface(s).norm;
    end
end

p1=p_start;
norm1=surface(s).norm;
p=surface(s).vert(1,:);
n_ref=0;
alpha=0;
Ra=1;
s2_t=0;

% plot3(p1(1),p1(2),p1(3),'ro','MarkerFaceColor','r');
% hold on


while alpha < Ra 
    n_ref=n_ref+1;

    [Rot] = rot_mat_emiss(norm1,p1,p);
    [omega] = random_dir(Rot);
    d_old=1000000;

    

    for m=1:1:length(surface(s).match)
        %
        s2=surface(s).match(m);

        vert_s=surface(s2).vert;
        n_plane=surface(s2).norm;
        p_plane=surface(s2).center;
        area_s=surface(s2).area;  
        %
        [p_end,d] = line_plane_inters(p1,omega,n_plane,p_plane);
        [flag1] = is_inside2(area_s,vert_s(:,1),vert_s(:,2),vert_s(:,3),p_end);
        if flag1==1 && d>0 && abs(d)<abs(d_old)
            
            d_old=d;
            s2_t=s2;
            p_rec=p_end;
            
        end
        
    end
    if exist("p_rec")==1   
        
        ref(n_ref).p=[p1;p_rec];
        s=s2_t;
        if strcmp(surface(s).item,'board')==1
           retta=(p1-p_rec)/norm(p1-p_rec);
            if dot(retta,surface(s).norm) < 0
                norm1=-surface(s).norm;
            else 
                norm1=surface(s).norm;
            end
        else
            norm1=surface(s).norm;
        end
%         plot3([p1(1);p_rec(1)],[p1(2);p_rec(2)],[p1(3);p_rec(3)],'y','LineWidth',5);
        p1=p_rec;
        
        p=surface(s).vert(1,:);     
        alpha=surface(s).prop_opt(2);
        Ra=rand;
        clear p_rec
       
    else
        disp('NO')
        flag3=1;
        ref=[];
        break
%         Om=p1+omega'*20;
%         plot3([p1(:,1);Om(:,1)],[p1(:,2);Om(:,2)],[p1(:,3);Om(:,3)],'linewidth',5)
    %% 
    %% 
    end
end

s_end=s2_t;
p_end=p1;
n_end=n_ref;
end