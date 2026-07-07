function [elem,total_nodes,Con] = node_face_creator3(C,L,n,An,th,item)

r1=[1 0 0; 0 cosd(An(1)) -sind(An(1)); 0 sind(An(1)) cosd(An(1))];
r2=[cosd(An(2)) 0 sind(An(2)); 0 1 0; -sind(An(2)) 0 cosd(An(2))];
r3=[cosd(An(3)) -sind(An(3)) 0; sind(An(3)) cosd(An(3)) 0; 0 0 1];
rot=r3*r2*r1;
dx=(L(1)/(n(1)-1))/2; dy=(L(2)/(n(2)-1))/2;

total_nodes=n(1)*n(2);

x=linspace(-L(1)/2,L(1)/2,n(1));
y=linspace(-L(2)/2,L(2)/2,n(2));
[X,Y]=meshgrid(x,y);
Z=zeros(size(X));

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

points=points*rot;

a=dy*dx;

Con=zeros(total_nodes,total_nodes);

kk = @(j,i) j+(i-1)*n(1);

k=1;
for i=1:1:n(2)
    for j=1:1:n(1)   
         
%         elem(k).ID=k;
%         elem(k).item=item;
        elem(k).node=[X(i,j) Y(i,j) Z(i,j)]*rot+C;

        if i==1
            if j==1
                %case1
                elem(k).type='v';
                id=[1,8,9,2];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a;
                elem(k).Ac=[dy,dx,0,0]*th;
                elem(k).V=elem(k).Af*th;
                elem(k).node_diff=elem(k).node;
                Con(kk(j,i),kk(j+1,i))=1;
                Con(kk(j,i),kk(j,i+1))=2; 
                k=k+1;
          
            elseif j==n(1)
                %case2
                elem(k).type='v';
                id=[1,6,7,8];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a;
                elem(k).Ac=[0,dx,dy,0]*th;
                elem(k).V=elem(k).Af*th;                
                Con(kk(j,i),kk(j,i+1))=2;
                Con(kk(j,i),kk(j-1,i))=3;
                elem(k).node_diff=elem(k).node;
                k=k+1;
                
            else
                %case3
                elem(k).type='s';
                id=[6,7,9,2];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a*2;
                elem(k).Ac=[dy,dx*2,dy,0]*th;
                elem(k).V=elem(k).Af*th;
                Con(kk(j,i),kk(j+1,i))=1;
                Con(kk(j,i),kk(j,i+1))=2;
                Con(kk(j,i),kk(j-1,i))=3;
                elem(k).node_diff=elem(k).node;
                k=k+1;
                
            end
        elseif i==n(2)
            if j==1
                %case4
                elem(k).type='v';
                id=[1,2,3,4];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a;
                elem(k).Ac=[dy,0,0,dx]*th;
                elem(k).V=elem(k).Af*th; 
                Con(kk(j,i),kk(j+1,i))=1;                
                Con(kk(j,i),kk(j,i-1))=4;
                elem(k).node_diff=elem(k).node;
                k=k+1;

            elseif j==n(1)
                %case5
                elem(k).type='v';
                id=[1,4,5,6];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a;
                elem(k).Ac=[0,0,dy,dx]*th;
                elem(k).V=elem(k).Af*th;
                Con(kk(j,i),kk(j-1,i))=3;
                Con(kk(j,i),kk(j,i-1))=4;
                elem(k).node_diff=elem(k).node;
                k=k+1;
               
            else
                %case6
                elem(k).type='s';
                id=[2,3,5,6];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a*2;
                elem(k).Ac=[dy,0,dy,dx*2]*th;
                elem(k).V=elem(k).Af*th;
                Con(kk(j,i),kk(j+1,i))=1;               
                Con(kk(j,i),kk(j-1,i))=3;
                Con(kk(j,i),kk(j,i-1))=4;
                elem(k).node_diff=elem(k).node;
                k=k+1;
                
            end
        else
            if j==1
                %case7
                elem(k).type='s';
                id=[4,8,9,3];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a*2;
                elem(k).Ac=[dy*2,dx,0,dx]*th;
                elem(k).V=elem(k).Af*th;
                Con(kk(j,i),kk(j+1,i))=1;
                Con(kk(j,i),kk(j,i+1))=2;          
                Con(kk(j,i),kk(j,i-1))=4;
                elem(k).node_diff=elem(k).node;
                k=k+1;
                
            elseif j==n(1)
                %case8 
                elem(k).type='s';
                id=[4,5,7,8];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a*2;
                elem(k).Ac=[0,dx,dy*2,dx]*th;
                elem(k).V=elem(k).Af*th;
                Con(kk(j,i),kk(j,i+1))=2;
                Con(kk(j,i),kk(j-1,i))=3;
                Con(kk(j,i),kk(j,i-1))=4;
                elem(k).node_diff=elem(k).node;
                k=k+1;
                
            else
                %case9
                elem(k).type='c';
                id=[3,5,7,9];
                elem(k).vertf=points(id,:)+elem(k).node;
                elem(k).Af=a*4;
                elem(k).Ac=[dy,dx,dy,dx]*2*th;
                elem(k).V=elem(k).Af*th;
                Con(kk(j,i),kk(j+1,i))=1;
                Con(kk(j,i),kk(j,i+1))=2;
                Con(kk(j,i),kk(j-1,i))=3;
                Con(kk(j,i),kk(j,i-1))=4;
                elem(k).node_diff=elem(k).node;
                k=k+1;
            
            end
        end
    end
end

% Connect=sparse(Con);
% i=3;
% figure
%   for i=1:1:length(elem)
% 
% plot3(elem(i).node(1),elem(i).node(2),elem(i).node(3),'o','Color','b','MarkerSize',10,...
%     'MarkerFaceColor','b')
%  hold on
% patch(elem(i).vert(:,1),elem(i).vert(:,2),elem(i).vert(:,3),'y')
%  hold on
% 
%  end
% 
%  
% 
% figure
% spy(Connect)


end