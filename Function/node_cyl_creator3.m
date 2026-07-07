function [elem,Con] = node_cyl_creator3(Nodes3D,Central,Bricks,R,L,Nt,Nr,Nz,total_nodes)

nc=size(Central,1);
nb=size(Bricks,1);
n=nc+nb;

[Afb,Aft,Afl,Acb,Act,Acl,V] = cylinder_areas(R,Nr,Nt,L,Nz);
Afbc = Afb(1,:)*Nt;
Aftc = Aft(1,:)*Nt;
Acbc = Acb(1,:);
% bottom, central
Acbc(3) = Acbc(3)*Nt;
Acbc(1) = 0;
Acbc(4) = 0;
% top, central
Actc = Act(1,:);
Actc(6) = Actc(6)*Nt;
Actc(1) = 0;
Actc(4) = 0;
% internal, central
Aclc = Acl(1,:);
Aclc(1) = 0;
Aclc(4) = 0;
Aclc(3) = Aclc(3)*Nt;
Aclc(6) = Aclc(6)*Nt;

elem=struct('ID',cell(1,n),...
            'ex_in',cell(1,n),...
            'item',cell(1,n),...
            'number',cell(1,n),...
            'ID_item',cell(1,n),...
            'node',cell(1,n),...
            'node_diff',cell(1,n),...
            'type',cell(1,n),...
            'face',cell(1,n),...
            'vertf',cell(1,n), ...
            'Af',cell(1,n), ...
            'Ac',cell(1,n), ...
            'V',cell(1,n),...
            'prop_mech',cell(1,n));
            %'prop_opt',cell(1,n));


Con= zeros(total_nodes,total_nodes);

% k=@(i,j,h) i+(j-1)*Nt+(h-1)*Nt*Nr;
k=@(i,j,h) (Nt*(Nr-1) + 1)*(h-1) + (j-2)*Nt*(j>1) + 1 + i*(j>1);

for h=1:1:(Nz-1) 
    for j=1:1:Nr
        if j == 1
            i = 1;
            m=k(i,j,h);                       
            if h == 1
                 elem(m).face=2;
                 elem(m).type='cb';
                 elem(m).node=mean(Nodes3D(Central(h,1:Nt),:),1);
                 elem(m).vertf=Nodes3D(Central(h,1:Nt),:);
                 elem(m).node_diff=mean(Nodes3D(Central(h,:)),1);
                 %case1
                 Con(m,k(1:Nt,j+1,h))=2; % Nt connection (+j dir)
                 Con(m,k(i,j,h+1))=3;

                 elem(m).Af=Afbc;
                 elem(m).Ac=Acbc;
                 elem(m).V=V(j)*Nt;

            elseif h == (Nz-1)
                 elem(m).face=1;
                 elem(m).type='ct';
                 elem(m).node=mean(Nodes3D(Central(h,(Nt+1):end),:),1);
                 elem(m).vertf=Nodes3D(Central(h,(Nt+1):end),:);
                 elem(m).node_diff=mean(Nodes3D(Central(h,:)),1);
                 %case2
                 Con(m,k(1:Nt,j+1,h))=2; % Nt connection (+j dir)
                 Con(m,k(i,j,h-1))=6;

                 elem(m).Af=Aftc;
                 elem(m).Ac=Actc;
                 elem(m).V=V(j)*Nt;

            else
                 elem(m).type='i';
                 elem(m).node=mean(Nodes3D(Central(h,:),:),1);
                 elem(m).node_diff=mean(Nodes3D(Central(h,:)),1);
                 %case3
                 Con(m,k(1:Nt,j+1,h))=2; % Nt connection (+j dir)
                 Con(m,k(i,j,h-1))=6;
                 Con(m,k(i,j,h+1))=3;

                 elem(m).Af=zeros(1,6);
                 elem(m).Ac=Aclc;
                 elem(m).V=V(j)*Nt;
            end
        else
            for i=1:1:Nt
                m=k(i,j,h);
                
                if h==1
                 elem(m).face=2;
                 
                 
%                  if j==1
%                      p=i+(h-1)*Nt;
%                      elem(m).type='ct';
%                      elem(m).node=mean(Nodes3D(Prisms(p,1:3),:),1);
%                      elem(m).vertf(1:3,:)=Nodes3D(Prisms(p,1:3),:);
%                      
%                 
%                      if i==1  
%                          %case1                     
%                      Con(m,k(i+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;
%                      Con(m,k(i,j,h+1))=3;
%                      Con(m,k(Nt,j,h))=4;
%                 
%                      elseif i==Nt
%                          %case2
%                      Con(m,k(i-Nt+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;
%                      Con(m,k(i,j,h+1))=3;
%                      Con(m,k(i-1,j,h))=4;
%                 
%                      else
%                          %case3
%                      Con(m,k(i+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;
%                      Con(m,k(i,j,h+1))=3;
%                      Con(m,k(i-1,j,h))=4;                   
%                      end
%                 
                 if j==Nr
                     elem(m).face=[2,i+2];
                     elem(m).type='s';
                     b=i+(j-2)*Nt+(h-1)*Nt*(Nr-1);
                     elem(m).node=mean(Nodes3D(Bricks(b,1:2),:),1);
                     elem(m).vertf(1:4,:)=Nodes3D(Bricks(b,1:4),:);
                     elem(m).vertf(5,:)=Nodes3D(Bricks(b,1),:);
                     elem(m).vertf(6,:)=Nodes3D(Bricks(b,2),:);
                     elem(m).vertf(7,:)=Nodes3D(Bricks(b,6),:);
                     elem(m).vertf(8,:)=Nodes3D(Bricks(b,5),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);
                
                     if i==1
                         %case4
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(Nt,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                
                     elseif i==Nt
                         %case5
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                         
                     else
                         %case6
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     end
                
                 else                   
                     b=i+(j-2)*Nt+(h-1)*Nt*(Nr-1);
                     elem(m).type='cq';
                     elem(m).node=mean(Nodes3D(Bricks(b,1:4),:),1);
                     elem(m).vertf(1:4,:)=Nodes3D(Bricks(b,1:4),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);
                
                     if i==1
                         %case7
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(Nt,j,h))=4; 
                     Con(m,k(i,j-1,h))=5;
                
                     elseif i==Nt
                         %case8
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4; 
                     Con(m,k(i,j-1,h))=5;
                         
                     else
                         %case9
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4; 
                     Con(m,k(i,j-1,h))=5;
                     end
                 end
                 elem(m).Af=Afb(j,:);
                 elem(m).Ac=Acb(j,:);
                 elem(m).V=V(j);
                
                elseif h==Nz-1
                  elem(m).face=1;
%                  if j==1
%                      p=i+(h-1)*Nt;
%                      elem(m).type='ct';
%                      elem(m).node=mean(Nodes3D(Prisms(p,4:6),:),1);
%                      elem(m).vertf(1:3,:)=Nodes3D(Prisms(p,4:6),:); 
%                 
%                      if i==1
%                          %case10
%                      Con(m,k(i+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;                     
%                      Con(m,k(Nt,j,h))=4;
%                      Con(m,k(i,j,h-1))=6;
%                 
%                      elseif i==Nt
%                          %case11
%                      Con(m,k(i-Nt+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;                     
%                      Con(m,k(i-1,j,h))=4;
%                      Con(m,k(i,j,h-1))=6;
%                          
%                      else
%                          %case12
%                      Con(m,k(i+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;                     
%                      Con(m,k(i-1,j,h))=4;
%                      Con(m,k(i,j,h-1))=6;
%                      end
                 if j==Nr
                     elem(m).face=[1,i+2];
                     elem(m).type='s';
                     b=i+(j-2)*Nt+(h-1)*Nt*(Nr-1);
                     elem(m).node=mean(Nodes3D(Bricks(b,5:6),:),1);
                     elem(m).vertf(1:4,:)=Nodes3D(Bricks(b,5:8),:);
                     elem(m).vertf(5,:)=Nodes3D(Bricks(b,5),:);
                     elem(m).vertf(6,:)=Nodes3D(Bricks(b,6),:);
                     elem(m).vertf(7,:)=Nodes3D(Bricks(b,2),:);
                     elem(m).vertf(8,:)=Nodes3D(Bricks(b,1),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);
                
                     if i==1
                         %case13
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(Nt,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                
                     elseif i==Nt
                         %case14
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                         
                     else
                         %case15
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                     end
                 else
                     b=i+(j-2)*Nt+(h-1)*Nt*(Nr-1);
                     elem(m).type='cq';
                     elem(m).node=mean(Nodes3D(Bricks(b,5:8),:),1);
                     elem(m).vertf(1:4,:)=Nodes3D(Bricks(b,5:8),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);
                
                     if i==1
                         %case16
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;                     
                     Con(m,k(Nt,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                
                     elseif i==Nt
                         %case17
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;                     
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                         
                     else
                         %case18
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;                     
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                     end
                 end
                 elem(m).Af=Aft(j,:);
                 elem(m).Ac=Act(j,:);
                 elem(m).V=V(j);
                
                else
                 elem(m).face=[];
%                  if j==1
%                      p=i+(h-1)*Nt;
%                      elem(m).type='i';
%                      elem(m).node=mean(Nodes3D(Prisms(p,:),:),1);                     
%                 
%                      if i==1
%                          %case19
%                      Con(m,k(i+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;
%                      Con(m,k(i,j,h+1))=3;
%                      Con(m,k(Nt,j,h))=4;
%                      Con(m,k(i,j,h-1))=6;
%                 
%                      elseif i==Nt
%                          %case20
%                      Con(m,k(i-Nt+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;
%                      Con(m,k(i,j,h+1))=3;
%                      Con(m,k(i-1,j,h))=4;
%                      Con(m,k(i,j,h-1))=6;
%                          
%                      else
%                          %case21
%                      Con(m,k(i+1,j,h))=1;
%                      Con(m,k(i,j+1,h))=2;
%                      Con(m,k(i,j,h+1))=3;
%                      Con(m,k(i-1,j,h))=4;
%                      Con(m,k(i,j,h-1))=6;
%                      end
%                 
                 if j==Nr
                     elem(m).face=i+2;
                     elem(m).type='cq';
                     b=i+(j-2)*Nt+(h-1)*Nt*(Nr-1);
                     elem(m).node=mean(Nodes3D(Bricks(b,[1,2,5,6]),:),1);
                     elem(m).vertf(1,:)=Nodes3D(Bricks(b,1),:);
                     elem(m).vertf(2,:)=Nodes3D(Bricks(b,2),:);
                     elem(m).vertf(3,:)=Nodes3D(Bricks(b,6),:);
                     elem(m).vertf(4,:)=Nodes3D(Bricks(b,5),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);
                
                     if i==1
                         %case22
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(Nt,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                
                     elseif i==Nt
                         %case23
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                         
                     else
                         %case24
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                     end
                
                 else
                     b=i+(j-2)*Nt+(h-1)*Nt*(Nr-1);
                     elem(m).type='i';
                     elem(m).node=mean(Nodes3D(Bricks(b,:),:),1); 
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);
                
                     if i==1
                         %case25
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(Nt,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                
                     elseif i==Nt
                         %case26
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                         
                     else
                         %case27
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     Con(m,k(i,j-1,h))=5;
                     Con(m,k(i,j,h-1))=6;
                     end
                
                 end
                 elem(m).Af=Afl(j,:);
                 elem(m).Ac=Acl(j,:);
                 elem(m).V=V(j);
                end
            end
        
        end
    end
end

end

