function [elem,total_nodes,Con] = node_box_creator2(C,L,n,An,th)

th1=th(1); th2=th(2); th3=th(3); th4=th(4); th5=th(5); th6=th(6);
r1=[1 0 0; 0 cosd(An(1)) -sind(An(1)); 0 sind(An(1)) cosd(An(1))];
r2=[cosd(An(2)) 0 sind(An(2)); 0 1 0; -sind(An(2)) 0 cosd(An(2))];
r3=[cosd(An(3)) -sind(An(3)) 0; sind(An(3)) cosd(An(3)) 0; 0 0 1];
rot=r3*r2*r1;
dx=(L(1)/(n(1)-1))/2; dy=(L(2)/(n(2)-1))/2; dz=(L(3)/(n(3)-1))/2;

total_nodes=(4+((n(1)-2)*2)+((n(2)-2)*2)+((n(1)-2)*(n(2)-2)))*2+...
              (n(3)-2)*4+...
              ((n(3)-2)*(n(1)-2))*2+...
              ((n(3)-2)*(n(2)-2))*2;

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
            'prop_mech',cell(1,total_nodes));
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

points=points*rot;
% figure
% plot3(points(:,1),points(:,2),points(:,3),'o')

Con=zeros(total_nodes,total_nodes);
A1=@(j,i,h) (i-1)*n(1);
B1=@(j,i,h) n(1)*n(2) + ((n(1)-1)*2 + (n(2)-1)*2)*(h-2);
C1=@(j,i,h) 2*(i-2) + n(1);
D1=@(j,i,h) j;
E1=@(j,i,h) 2;

kk=@(j,i,h) A1(j,i,h)*(h==1 || h==n(3)) +...
    B1(j,i,h)*(h~=1) +...
    C1(j,i,h)*(h~=1 && h~=n(3) && i~=1) +...
    D1(j,i,h)*(~(h~=1 && h~=n(3) && i~=1 && i~=n(2) && j==n(1))) +...
    E1(j,i,h)*(h~=1 && h~=n(3) && i~=1 && i~=n(2) && j==n(1));

k=1;
e=1;
for h=1:1:n(3)
    for i=1:1:n(2)
        for j=1:1:n(1)
%             elem(k).ID=k;
%             elem(k).item='box';
            elem(k).node=[X(i,j,h) Y(i,j,h) Z(i,j,h)]*rot+C;

            if h==1
                if i==1
                    if j==1
                        %case1
                        elem(k).face=[4,1,5];
                        elem(k).type='v';
                        id=[1,10,17,8,1,2,11,10,1,8,9,2];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz,dx*dz,dx*dy,dy*dz,dx*dz,dx*dy];
                        elem(k).Ac=[dz*th1+dy*th5,dz*th4+dx*th5,dx*th1+dy*th4,0,0,0];
                        elem(k).V=dx*th1*dz+dy*th4*dz+dx*dy*th5;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                        %case2
                        elem(k).face=[2,1,5];
                        elem(k).type='v';
                        id=[1,10,17,8,1,10,15,6,1,6,7,8];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz,dx*dz,dx*dy,dy*dz,dx*dz,dx*dy];
                        elem(k).Ac=[0,dz*th2+dx*th5,dx*th1+dy*th2,dz*th1+dx*th5,0,0];
                        elem(k).V=dx*th1*dz+dy*th2*dz+dx*dy*th5;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    else
                        %case3
                        elem(k).face=[1,5];
                        elem(k).type='s';
                        id=[2,11,15,6,2,6,7,9];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[0,dx*2*dz,dx*2*dy,0,dx*2*dz,dx*2*dy];
                        elem(k).Ac=[dz*th1+dy*th5,dx*2*th5,dx*2*th1,dz*th1+dy*th5,0,0];
                        elem(k).V=dx*2*th1*dz+dx*2*th5*dy;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    end
                elseif i==n(2)
                    if j==1
                        %case4
                        elem(k).face=[4,3,5];
                        elem(k).type='v';
                        id=[1,4,13,10,1,10,11,2,1,2,3,4];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz,dx*dz,dx*dy,dy*dz,dx*dz,dx*dy];
                        elem(k).Ac=[dz*th3+dx*th5,0,dy*th4+dx*th3,0,dz*th4+dx*th5,0];
                        elem(k).V=dx*th3*dz+dy*th4*dz+dx*dy*th5;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                        %case5
                        elem(k).face=[2,3,5];
                        elem(k).type='v';
                        id=[1,4,13,10,1,6,15,10,1,4,5,6];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz,dx*dz,dx*dy,dy*dz,dx*dz,dx*dy];
                        elem(k).Ac=[0,0,dy*th2+dx*th3,dz*th3+dy*th5,dz*th2+dx*th5,0];
                        elem(k).V=dx*th3*dz+dy*th2*dz+dx*dy*th5;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    else
                        %case6
                        elem(k).face=[3,5];
                        elem(k).type='s';
                        id=[2,11,15,6,2,3,5,6];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[0,dx*2*dz,dx*2*dy,0,dx*2*dz,dx*2*dy];
                        elem(k).Ac=[dz*th3+dy*th5,0,dx*2*th3,dz*th3+dy*th5,dx*2*th5,0];
                        elem(k).V=dx*2*th3*dz+dx*2*th5*dy;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    end
                else
                    if j==1
                        %case7
                        elem(k).face=[4,5];
                        elem(k).type='s';
                        id=[4,13,17,8,4,8,9,3];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*2*dz,0,dy*2*dx,dy*2*dz,0,dy*2*dx];
                        elem(k).Ac=[dy*2*th5,dz*th4+dx*th5,dy*2*th4,0,dz*th4+dx*th5,0];
                        elem(k).V=dy*2*dx*th5+dy*2*dz*th4;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                        %case8
                        elem(k).face=[2,5];
                        elem(k).type='s';
                        id=[4,8,17,13,4,8,7,5];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*2*dz,0,dy*2*dx,dy*2*dz,0,dy*2*dx];
                        elem(k).Ac=[0,dz*th2+dx*th5,dy*2*th2,dy*2*th5,dz*th2+dx*th5,0];
                        elem(k).V=dy*2*dx*th5+dy*2*dz*th2;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    else
                        %case9
                        elem(k).face=5;
                        elem(k).type='c';
                        id=[3,5,7,9];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[0,0,dx*2*dy*2,0,0,dx*2*dy*2];
                        elem(k).Ac=[dy*2*th5,dx*2*th5,0,dy*2*th5,dx*2*th5,0];
                        elem(k).V=dy*2*dx*2*th5;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                   
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    end
                end
            elseif h==n(3)
                if i==1
                    if j==1
                        %case10
                        elem(k).face=[4,1,6];
                        elem(k).type='v';
                        id=[1,8,26,19,1,2,20,19,1,8,9,2];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz,dx*dz,dx*dy,dy*dz,dx*dz,dx*dy];
                        elem(k).Ac=[dz*th1+dy*th6,dz*th4+dx*th6,0,0,0,dx*th1+dy*th4];
                        elem(k).V=dx*th1*dz+dy*th4*dz+dx*dy*th6;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                       
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                        %case11
                        elem(k).face=[2,1,6];
                        elem(k).type='v';
                        id=[1,8,26,19,1,6,24,19,1,8,7,6];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz,dx*dz,dx*dy,dy*dz,dx*dz,dx*dy];
                        elem(k).Ac=[0,dz*th2+dx*th6,0,dz*th1+dx*th6,0,dx*th1+dy*th2];
                        elem(k).V=dx*th1*dz+dy*th2*dz+dx*dy*th6;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                    
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                   
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    else
                        %case12
                        elem(k).face=[1,6];
                        elem(k).type='s';
                        id=[2,6,24,20,2,6,7,9];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[0,dx*2*dz,dx*2*dy,0,dx*2*dz,dx*2*dy];
                        elem(k).Ac=[dz*th1+dy*th6,dx*2*th6,0,dz*th1+dy*th6,0,dx*2*th1];
                        elem(k).V=dx*2*th1*dz+dx*2*th6*dy;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                        
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                        
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    end
                elseif i==n(2)
                    if j==1
                        %case13
                        elem(k).face=[4,3,6];
                        elem(k).type='v';
                        id=[1,4,22,19,1,2,20,19,1,4,3,2];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz,dx*dz,dx*dy,dy*dz,dx*dz,dx*dy];
                        elem(k).Ac=[dz*th3+dy*th6,0,0,0,dz*th4+dx*th6,dx*th3+dy*th4];
                        elem(k).V=dx*th3*dz+dy*th4*dz+dx*dy*th6;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                    
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6; 
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                        %case14
                        elem(k).face=[2,3,6];
                        elem(k).type='v';
                        id=[1,4,22,19,1,6,24,19,1,4,5,6];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz,dx*dz,dx*dy,dy*dz,dx*dz,dx*dy];
                        elem(k).Ac=[0,0,0,dz*th3+dy*th6,dz*th2+dx*th6,dx*th3+dy*th2];
                        elem(k).V=dx*th3*dz+dy*th2*dz+dx*dy*th6;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    else
                        %case15
                        elem(k).face=[3,6];
                        elem(k).type='s';
                        id=[2,6,24,20,2,6,5,3];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[0,dx*2*dz,dx*2*dy,0,dx*2*dz,dx*2*dy];
                        elem(k).Ac=[dz*th3+dy*th6,0,0,dz*th3+dy*th6,dx*2*th6,dx*2*th3];
                        elem(k).V=dx*2*th3*dz+dx*2*th6*dy;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                        
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    end
                else
                    if j==1
                        %case16
                        elem(k).face=[4,6];
                        elem(k).type='s';
                        id=[4,22,26,8,4,8,9,3];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*2*dz,0,dy*2*dx,dy*2*dz,0,dy*2*dx];
                        elem(k).Ac=[dy*2*th6,dz*th4+dx*th6,0,0,dz*th4+dx*th6,dy*2*th4];
                        elem(k).V=dy*2*dx*th6+dy*2*dz*th4;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                    
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                        %case17
                        elem(k).face=[2,6];
                        elem(k).type='s';
                        id=[4,8,26,22,4,8,7,5];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*2*dz,0,dy*2*dx,dy*2*dz,0,dy*2*dx];
                        elem(k).Ac=[0,dz*th2+dx*th6,0,dy*2*th6,dz*th2+dx*th6,dy*2*th2];
                        elem(k).V=dy*2*dx*th6+dy*2*dz*th2;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                   
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    else
                        %case18
                        elem(k).face=6;
                        elem(k).type='c';
                        id=[3,5,7,9];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[0,0,dx*2*dy*2,0,0,dx*2*dy*2];
                        elem(k).Ac=[dy*2*th6,dx*2*th6,0,dy*2*th6,dx*2*th6,0];
                        elem(k).V=dy*2*dx*2*th6;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;                    
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    end
                end
            else
                if i==1
                    if j==1
                        %case19
                        elem(k).face=[4,1];
                        elem(k).type='s';
                        id=[10,17,26,19,10,11,20,19];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz*2,dx*dz*2,0,dy*dz*2,dx*dz*2,0];
                        elem(k).Ac=[dz*2*th1,dz*2*th4,dx*th1+dy*th4,0,0,dx*th1+dy*th4];
                        elem(k).V=dz*2*th4*dy+dz*2*th1*dx;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;                        
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                        %case20
                        elem(k).face=[2,1];
                        elem(k).type='s';
                        id=[10,17,26,19,10,15,24,19];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz*2,dx*dz*2,0,dy*dz*2,dx*dz*2,0];
                        elem(k).Ac=[0,dz*2*th2,dx*th1+dy*th2,dz*2*th1,0,dx*th1+dy*th2];
                        elem(k).V=dz*2*th4*dy+dz*2*th1*dx;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                       
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    else
                        %case21
                        elem(k).face=1;
                        elem(k).type='c';
                        id=[11,15,24,20];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[0,dx*2*dz*2,0,0,dx*2*dz*2,0];
                        elem(k).Ac=[dz*2*th1,0,dx*2*th1,dz*th1*2,0,dx*th1*2];
                        elem(k).V=dz*2*dx*2*th1;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                       
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                        
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    end
                
                elseif i==n(2)
                    if j==1
                        %case22
                        elem(k).face=[4,3];
                        elem(k).type='s';
                        id=[10,13,22,19,10,11,20,19];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz*2,dx*dz*2,0,dy*dz*2,dx*dz*2,0];
                        elem(k).Ac=[dz*2*th3,0,dx*th3+dy*th4,0,dz*2*th4,dx*th3+dy*th4];
                        elem(k).V=dz*2*th4*dy+dz*2*th3*dx;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                       
                        Con(kk(j,i,h),kk(j,i,h+1))=3;                       
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                        %case23
                        elem(k).face=[2,3];
                        elem(k).type='s';
                        id=[10,13,22,19,10,15,24,19];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*dz*2,dx*dz*2,0,dy*dz*2,dx*dz*2,0];
                        elem(k).Ac=[0,0,dx*th3+dy*th2,dz*2*th3,dz*2*th2,dx*th3+dy*th2];
                        elem(k).V=dz*2*th2*dy+dz*2*th3*dx;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    else
                        %case24
                        elem(k).face=3;
                        elem(k).type='c';
                        id=[11,15,24,20];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[0,dx*2*dz*2,0,0,dx*2*dz*2,0];
                        elem(k).Ac=[dz*2*th3,0,dx*2*th3,dz*th3*2,0,dx*th3*2];
                        elem(k).V=dz*2*dx*2*th3;
                        Con(kk(j,i,h),kk(j+1,i,h))=1;                     
                        Con(kk(j,i,h),kk(j,i,h+1))=3;
                        Con(kk(j,i,h),kk(j-1,i,h))=4;                        
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    end
                else
                    if j==1
                     %case25
                        elem(k).face=4;
                        elem(k).type='c';
                        id=[13,17,26,22];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*2*dz*2,0,0,dy*2*dz*2,0,0];
                        elem(k).Ac=[0,dz*2*th4,dy*2*th4,0,dz*2*th4,dy*th4*2];
                        elem(k).V=dz*2*dy*2*th4;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;                
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;

                    elseif j==n(1)
                     %case26
                        elem(k).face=2;
                        elem(k).type='c';
                        id=[13,17,26,22];
                        elem(k).vertf=points(id,:)+elem(k).node;
                        elem(k).Af=[dy*2*dz*2,0,0,dy*2*dz*2,0,0];
                        elem(k).Ac=[0,dz*2*th2,dy*2*th2,0,dz*2*th2,dy*th2*2];
                        elem(k).V=dz*2*dy*2*th2;
                        Con(kk(j,i,h),kk(j,i+1,h))=2;
                        Con(kk(j,i,h),kk(j,i,h+1))=3;                   
                        Con(kk(j,i,h),kk(j,i-1,h))=5;
                        Con(kk(j,i,h),kk(j,i,h-1))=6;
                        elem(k).node_diff=elem(k).node;
                        k=k+1;
                    else
                    % case27
                       elim(e)=kk(i,j,h);
                       e=e+1;

                    end
               end
            end
        end
    end
end

% Con(elim,:)=[];
% Con(:,elim)=[];



% Connect=sparse(Con);

%   i=26;  
%   figure
%   for i=1:1:length(elem)
%    plot3(elem(i).node(1),elem(i).node(2),elem(i).node(3),'o','Color','b','MarkerSize',10,...
%     'MarkerFaceColor','b')
%     hold on
%     for j=1:length(elem(i).face)
%     patch(elem(i).vertf((1:4)+(j-1)*4,1),elem(i).vertf((1:4)+(j-1)*4,2),elem(i).vertf((1:4)+(j-1)*4,3),'y')
%     hold on
%     end
%   end
% figure
%   spy(Connect)

  end

