function [elem,total_nodes,Con] = node_solid_creator2(C,L,n,An)


r1=[1 0 0; 0 cosd(An(1)) -sind(An(1)); 0 sind(An(1)) cosd(An(1))];
r2=[cosd(An(2)) 0 sind(An(2)); 0 1 0; -sind(An(2)) 0 cosd(An(2))];
r3=[cosd(An(3)) -sind(An(3)) 0; sind(An(3)) cosd(An(3)) 0; 0 0 1];
rot=r3*r2*r1;
dx=(L(1)/(n(1)-1))/2; dy=(L(2)/(n(2)-1))/2; dz=(L(3)/(n(3)-1))/2;

total_nodes=n(1)*n(2)*n(3);

x=linspace(-L(1)/2,L(1)/2,n(1));
y=linspace(-L(2)/2,L(2)/2,n(2));
z=linspace(-L(3)/2,L(3)/2,n(3));
[X, Y, Z]=meshgrid(x,y,z);

elem=struct('ID',cell(1,total_nodes),...
            'ex_in',cell(1,total_nodes),...
            'item',cell(1,total_nodes),...
            'number',cell(1,total_nodes),...
            'ID_item',cell(1,total_nodes),...
            'node',cell(1,total_nodes),...
            'node_diff',cell(1,total_nodes),...
            'type',cell(1,total_nodes),...
            'face',cell(1,total_nodes),...
            'vertf',cell(1,total_nodes),...            
            'Af',cell(1,total_nodes),...
            'Ac',cell(1,total_nodes),...
            'V',cell(1,total_nodes),...
            'prop_mech',cell(1,total_nodes),...
            'dz_local',cell(1,total_nodes));
            %'prop_opt',cell(1,total_nodes));

points(1,:)=[0,0,0];
points(2,:)=[dx,0,0];
points(3,:)=[+dx,-dy,0];
points(4,:)=[0,-dy,0];
points(5,:)=[-dx,-dy,0];
points(6,:)=[-dx,0,0];
points(7,:)=[-dx,+dy,0];
points(8,:)=[0,dy,0];
points(9,:)=[+dx,dy,0];
points(10:18,:)=points(1:9,:)+[0,0,dz];
points(19:27,:)=points(1:9,:)+[0,0,-dz];

ax=dy*dz;
ay=dx*dz;
az=dx*dy;

v=dx*dy*dz;

points=points*rot;
% figure
% plot3(points(:,1),points(:,2),points(:,3),'o')

Con=zeros(total_nodes,total_nodes);
kk=@(j,i,h) j+(i-1)*n(1)+(h-1)*n(1)*n(2);


k=1;
for h=1:1:n(3)
    for i=1:1:n(2)
        for j=1:1:n(1)
%             elem(k).ID=k;
%             elem(k).item='paral';
            elem(k).node=[X(i,j,h) Y(i,j,h) Z(i,j,h)]*rot+C;

            if h==1
                if i==1
                    if j==1
                        %case1
                        elem(k).face=[4,1,5];
                        elem(k).type='v';
                        idf=[1,10,17,8,1,2,11,10,1,2,9,8];                        
                        elem(k).vertf=points(idf,:)+elem(k).node;                       
                        elem(k).Af=[0,0,0,ax,ay,az];
                        elem(k).Ac=[ax,ay,az,0,0,0];
                        elem(k).V=v;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        idc=[1,2,9,8,10,11,16,17];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    elseif j==n(1)
                        %case2
                        elem(k).face=[2,1,5];
                        elem(k).type='v';
                        idf=[1,10,17,8,1,10,15,6,1,8,7,6];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax,0,0,0,ay,az];
                        elem(k).Ac=[0,ay,az,ax,0,0];
                        elem(k).V=v;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        idc=[1,6,7,8,10,15,16,17];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    else
                        %case3
                        elem(k).face=[1,5];
                        elem(k).type='s';
                        idf=[2,11,15,6,2,6,7,9];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,0,0,ay,az]*2;
                        elem(k).Ac=[ax,ay*2,az*2,ax,0,0];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        idc=[2,6,7,8,10,15,16,18];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    end
                elseif i==n(2)
                    if j==1
                        %case4
                        elem(k).face=[4,3,5];
                        elem(k).type='v';
                        idf=[1,4,13,10,1,2,11,10,1,4,3,2];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,ay,0,ax,0,az];
                        elem(k).Ac=[ax,0,az,0,ay,0];
                        elem(k).V=v;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        idc=[1,2,3,4,10,11,12,13];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    elseif j==n(1)
                        %case5
                        elem(k).face=[2,3,5];
                        elem(k).type='v';
                        idf=[1,4,13,10,1,6,15,10,1,4,5,6];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax,ay,0,0,0,az];
                        elem(k).Ac=[0,0,az,ax,ay,0];
                        elem(k).V=v;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        idc=[1,4,5,6,10,13,14,15];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    else
                        %case6
                        elem(k).face=[3,5];
                        elem(k).type='s';
                        idf=[2,11,15,6,2,3,5,6];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,ay,0,0,0,az]*2;
                        elem(k).Ac=[ax,0,az*2,ax,ay*2,0];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5; 
                        idc=[2,3,5,6,11,12,14,15];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    end
                else
                    if j==1
                        %case7
                        elem(k).face=[4,5];
                        elem(k).type='s';
                        idf=[4,13,17,8,4,8,9,3];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,0,ax,0,az]*2;
                        elem(k).Ac=[ax*2,ay,az*2,0,ay,0];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        idc=[3,4,8,9,12,13,17,18];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    elseif j==n(1)
                        %case8
                        elem(k).face=[2,5];
                        elem(k).type='s';
                        idf=[4,8,17,13,4,8,7,5];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax,0,0,0,0,az]*2;
                        elem(k).Ac=[0,ay,az*2,ax*2,ay,0];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        idc=[4,5,7,8,13,14,16,17];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    else
                        %case9
                        elem(k).face=5;
                        elem(k).type='c';
                        idf=[3,5,7,9];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,0,0,0,az*4];
                        elem(k).Ac=[ax,ay,az*2,ax,ay,0]*2;
                        elem(k).V=v*4;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        idc=[3,5,7,9,12,14,16,18];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    end
                end
            elseif h==n(3)
                if i==1
                    if j==1
                        %case10
                        elem(k).face=[4,1,6];
                        elem(k).type='v';
                        idf=[1,8,26,19,1,2,20,19,1,8,9,2];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,az,ax,ay,0];
                        elem(k).Ac=[ax,ay,0,0,0,az];
                        elem(k).V=v;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                       
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[1,2,9,8,20,19,26,27];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    elseif j==n(1)
                        %case11
                        elem(k).face=[2,1,6];
                        elem(k).type='v';
                        idf=[1,8,26,19,1,6,24,19,1,8,7,6];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax,0,az,0,ay,0];
                        elem(k).Ac=[0,ay,0,ax,0,az];
                        elem(k).V=v;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                    
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                   
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[1,6,7,8,19,24,25,26];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    else
                        %case12
                        elem(k).face=[1,6];
                        elem(k).type='s';
                        idf=[2,6,24,20,2,6,7,9];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,az,0,ay,0]*2;
                        elem(k).Ac=[ax,0,az*2,ax,ay*2,0];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                        
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                        
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[2,6,7,9,20,24,25,27];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    end
                elseif i==n(2)
                    if j==1
                        %case13
                        elem(k).face=[4,3,6];
                        elem(k).type='v';
                        idf=[1,4,22,19,1,2,20,19,1,4,3,2];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,ay,az,ax,0,0];
                        elem(k).Ac=[ax,0,0,0,ay,az];
                        elem(k).V=v;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                    
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[1,2,3,4,19,20,21,22];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    elseif j==n(1)
                        %case14
                        elem(k).face=[2,3,6];
                        elem(k).type='v';
                        idf=[1,4,22,19,1,6,24,19,1,4,5,6];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax,ay,az,0,0,0];
                        elem(k).Ac=[0,0,0,ax,ay,az];
                        elem(k).V=v;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[1,4,5,6,19,22,23,24];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    else
                        %case15
                        elem(k).face=[3,6];
                        elem(k).type='s';
                        idf=[2,6,24,20,2,6,5,3];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,ay,az,0,0,0]*2;
                        elem(k).Ac=[ax,0,0,ax,ay*2,az*2];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                        
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[2,3,5,6,20,21,23,24];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    end
                else
                    if j==1
                        %case16
                        elem(k).face=[4,6];
                        elem(k).type='s';
                        idf=[4,22,26,8,4,8,9,3];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,az,ax,0,0]*2;
                        elem(k).Ac=[ax*2,ay,0,0,ay,az*2];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                    
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[3,4,8,9,27,21,22,26];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    elseif j==n(1)
                        %case17
                        elem(k).face=[2,6];
                        elem(k).type='s';
                        idf=[4,8,26,22,4,8,7,5];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax,0,az,0,0,0]*2;
                        elem(k).Ac=[0,ay,0,ax*2,ay,az*2];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                   
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[4,5,7,8,22,23,25,26];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    else
                        %case18
                        elem(k).face=6;
                        elem(k).type='c';
                        idf=[3,5,7,9];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,az*4,0,0,0];
                        elem(k).Ac=[ax,ay,0,ax,ay,az*2]*2;
                        elem(k).V=v*4;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                    
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[3,5,7,9,21,23,25,27];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    end
                end
            else
                if i==1
                    if j==1
                        %case19
                        elem(k).face=[4,1];
                        elem(k).type='s';
                        idf=[10,17,26,19,10,11,20,19];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,0,ax,ay,0]*2;
                        elem(k).Ac=[ax*2,ay*2,az,0,0,az];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;                        
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[10,11,18,17,20,19,26,27];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    elseif j==n(1)
                        %case20
                        elem(k).face=[2,1];
                        elem(k).type='s';
                        idf=[10,17,26,19,10,15,24,19];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax,0,0,0,ay,0]*2;
                        elem(k).Ac=[0,ay*2,az,ax*2,0,az];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                       
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[10,15,16,17,19,24,25,26];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    else
                        %case21
                        elem(k).face=1;
                        elem(k).type='c';
                        idf=[11,15,24,20];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,0,0,ay*4,0];
                        elem(k).Ac=[ax,ay*2,az,ax,0,az]*2;
                        elem(k).V=v*4;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                        
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[10,15,16,18,20,24,25,27];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    end
                
                elseif i==n(2)
                    if j==1
                        %case22
                        elem(k).face=[4,3];
                        elem(k).type='s';
                        idf=[10,13,22,19,10,11,20,19];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,ay,0,ax,0,0]*2;
                        elem(k).Ac=[ax*2,0,az,0,ay*2,az];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                       
                        Con(kk(j,i,h),kk(j,i,h+1))=3;                       
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[11,10,13,12,20,21,22,19];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    elseif j==n(1)
                        %case23
                        elem(k).face=[2,3];
                        elem(k).type='s';
                        idf=[10,13,22,19,10,15,24,19];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax,ay,0,0,0,0]*2;
                        elem(k).Ac=[0,0,az,ax*2,ay*2,az];
                        elem(k).V=v*2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[10,13,14,15,19,22,23,24];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    else
                        %case24
                        elem(k).face=3;
                        elem(k).type='c';
                        idf=[11,15,24,20];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,ay*4,0,0,0,0];
                        elem(k).Ac=[ax,0,az,ax,ay*2,az]*2;
                        elem(k).V=v*4;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                     
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[11,12,14,15,20,21,23,24];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                    end
            else
                if j==1
                     %case25
                        elem(k).face=4;
                        elem(k).type='c';
                        idf=[13,17,26,22];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[0,0,0,ax*4,0,0];
                        elem(k).Ac=[ax*2,ay,az,0,ay,az]*2;
                        elem(k).V=v*4;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;                
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[13,14,16,17,22,26,27,21];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                elseif j==n(1)
                     %case26
                        elem(k).face=2;
                        elem(k).type='c';
                        idf=[13,17,26,22];
                        elem(k).vertf=points(idf,:)+elem(k).node;
                        elem(k).Af=[ax*4,0,0,0,0,0];
                        elem(k).Ac=[0,ay,az,ax*2,ay,az]*2;
                        elem(k).V=v*4;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        idc=[13,17,16,14,22,26,25,23];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;
                else
                     %case27
                        elem(k).face=[];
                        elem(k).type='i';                      
                        elem(k).vertf=[];
                        elem(k).Af=[0,0,0,0,0,0];
                        elem(k).Ac=[ax,ay,az,ax,ay,az]*4;
                        elem(k).V=v*8;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6; 
                        idc=[12,14,16,18,21,23,25,27];
                        elem(k).node_diff =mean(points(idc,:)+elem(k).node,1);
                        k=k+1;

                end
               end
            end
        end
    end
end

% Connect=sparse(Con);



% %   i=26;  
%   figure
%    for i=1:1:length(elem)
%    plot3(elem(i).node(1),elem(i).node(2),elem(i).node(3),'o','Color','b','MarkerSize',10,...
%     'MarkerFaceColor','b')
%     hold on
%     for j=1:length(elem(i).face)
%     patch(elem(i).vertf((1:4)+(j-1)*4,1),elem(i).vertf((1:4)+(j-1)*4,2),elem(i).vertf((1:4)+(j-1)*4,3),'y')
%     hold on
%     end
%    end
% figure
%    spy(Connect)
  end