function [elem,Con] = node_sphcap_creator(Nodes3D,Central,Bricks,R_int,R_out,h_cap,Ntheta,Nr_cap,Nt,total_nodes,phi_block_size)

if nargin<11 || isempty(phi_block_size)
    phi_block_size=1; % 1 = un gruppo di superficie per cella (esatto);
                       % Nt = un gruppo per l'intero anello. Vedi
                       % sphcap_face.m per come viene usato per generare
                       % le normali del guscio.
end
n_phi_blocks = ceil(Nt/phi_block_size);
phi_block = @(i) floor((i-1)/phi_block_size)+1;
% ID locali: l'apice (j==1) e' un unico nodo fuso -- un solo slot, NON
% n_phi_blocks (altrimenti gli altri restano vuoti, mai assegnati a
% nessun elemento -- crash a valle in surf_global.m). Anelli regolari
% (j=2..Ntheta) usano n_phi_blocks slot ciascuno, a partire da dopo lo
% slot dell'apice.
n_shell_slots = 1 + (Ntheta-1)*n_phi_blocks; % per guscio (interno o esterno)
face_inner = @(j,ii) (j==1)*1 + (j>1)*(1+(j-2)*n_phi_blocks+phi_block(ii));
face_outer = @(j,ii) n_shell_slots + (j==1)*1 + (j>1)*(1+(j-2)*n_phi_blocks+phi_block(ii));
face_lateral = @(ii) 2*n_shell_slots + ii;

nc=size(Central,1);
nb=size(Bricks,1);
n=nc+nb;

[Afb,Aft,Afl,Acb,Act,Acl,V,Geom_r] = sphcap_areas(R_int,R_out,h_cap,Ntheta,Nr_cap,Nt);

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
            'prop_mech',cell(1,n),...
            'dz_local',cell(1,n));

Con= zeros(total_nodes,total_nodes);

k=@(i,j,h) (Nt*(Ntheta-1) + 1)*(h-1) + (j-2)*Nt*(j>1) + 1 + i*(j>1);
bfun=@(i,j,h) i+(j-2)*Nt+(h-1)*Nt*(Ntheta-1);

for h=1:1:Nr_cap

    Afbc = reshape(Afb(1,h,:),1,6)*Nt;
    Aftc = reshape(Aft(1,h,:),1,6)*Nt;
    Acbc = reshape(Acb(1,h,:),1,6);
    Acbc(3) = Acbc(3)*Nt;
    Acbc(1) = 0;
    Acbc(4) = 0;
    Actc = reshape(Act(1,h,:),1,6);
    Actc(6) = Actc(6)*Nt;
    Actc(1) = 0;
    Actc(4) = 0;
    Aclc = reshape(Acl(1,h,:),1,6);
    Aclc(1) = 0;
    Aclc(4) = 0;
    Aclc(3) = Aclc(3)*Nt;
    Aclc(6) = Aclc(6)*Nt;

    for j=1:1:Ntheta
        if j == 1
            i = 1;
            m=k(i,j,h);
            if h == 1
                 elem(m).face=face_inner(j,1);
                 elem(m).type='cb';
                 elem(m).node=mean(Nodes3D(Central(h,1:Nt),:),1);
                 elem(m).vertf=Nodes3D(Central(h,1:Nt),:);
                 elem(m).node_diff=mean(Nodes3D(Central(h,:)),1);
                 Con(m,k(1:Nt,j+1,h))=2;
                 if h<Nr_cap
                 Con(m,k(i,j,h+1))=3;
                 end

                 elem(m).Af=Afbc;
                 elem(m).Ac=Acbc;
                 elem(m).V=V(j,h)*Nt;
                 elem(m).dz_local=Geom_r(j,h)*Nt;

            elseif h == Nr_cap
                 elem(m).face=face_outer(j,1);
                 elem(m).type='ct';
                 elem(m).node=mean(Nodes3D(Central(h,(Nt+1):end),:),1);
                 elem(m).vertf=Nodes3D(Central(h,(Nt+1):end),:);
                 elem(m).node_diff=mean(Nodes3D(Central(h,:)),1);
                 Con(m,k(1:Nt,j+1,h))=2;
                 Con(m,k(i,j,h-1))=6;

                 elem(m).Af=Aftc;
                 elem(m).Ac=Actc;
                 elem(m).V=V(j,h)*Nt;
                 elem(m).dz_local=Geom_r(j,h)*Nt;

            else
                 elem(m).type='i';
                 elem(m).node=mean(Nodes3D(Central(h,:),:),1);
                 elem(m).node_diff=mean(Nodes3D(Central(h,:)),1);
                 Con(m,k(1:Nt,j+1,h))=2;
                 Con(m,k(i,j,h-1))=6;
                 if h<Nr_cap
                 Con(m,k(i,j,h+1))=3;
                 end

                 elem(m).Af=zeros(1,6);
                 elem(m).Ac=Aclc;
                 elem(m).V=V(j,h)*Nt;
                 elem(m).dz_local=Geom_r(j,h)*Nt;
            end
        else
            for i=1:1:Nt
                m=k(i,j,h);

                if h==1
                 elem(m).face=face_inner(j,i);

                 if j==Ntheta
                     elem(m).face=[face_inner(j,i),face_lateral(i)];
                     elem(m).type='s';
                     b=bfun(i,j,h);
                     elem(m).node=mean(Nodes3D(Bricks(b,1:2),:),1);
                     elem(m).vertf(1:4,:)=Nodes3D(Bricks(b,1:4),:);
                     elem(m).vertf(5,:)=Nodes3D(Bricks(b,1),:);
                     elem(m).vertf(6,:)=Nodes3D(Bricks(b,2),:);
                     elem(m).vertf(7,:)=Nodes3D(Bricks(b,6),:);
                     elem(m).vertf(8,:)=Nodes3D(Bricks(b,5),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);

                     if i==1
                     Con(m,k(i+1,j,h))=1;
                     if h<Nr_cap
                     Con(m,k(i,j,h+1))=3;
                     end
                     Con(m,k(Nt,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end

                     elseif i==Nt
                     Con(m,k(i-Nt+1,j,h))=1;
                     if h<Nr_cap
                     Con(m,k(i,j,h+1))=3;
                     end
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end

                     else
                     Con(m,k(i+1,j,h))=1;
                     if h<Nr_cap
                     Con(m,k(i,j,h+1))=3;
                     end
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     end

                 else
                     b=bfun(i,j,h);
                     elem(m).type='cq';
                     elem(m).node=mean(Nodes3D(Bricks(b,1:4),:),1);
                     elem(m).vertf(1:4,:)=Nodes3D(Bricks(b,1:4),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);

                     if i==1
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     if h<Nr_cap
                     Con(m,k(i,j,h+1))=3;
                     end
                     Con(m,k(Nt,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end

                     elseif i==Nt
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     if h<Nr_cap
                     Con(m,k(i,j,h+1))=3;
                     end
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end

                     else
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     if h<Nr_cap
                     Con(m,k(i,j,h+1))=3;
                     end
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     end
                 end
                 elem(m).Af=reshape(Afb(j,h,:),1,6);
                 elem(m).Ac=reshape(Acb(j,h,:),1,6);
                 elem(m).V=V(j,h);
                 elem(m).dz_local=Geom_r(j,h);

                elseif h==Nr_cap
                  elem(m).face=face_outer(j,i);
                 if j==Ntheta
                     elem(m).face=[face_outer(j,i),face_lateral(i)];
                     elem(m).type='s';
                     b=bfun(i,j,h);
                     elem(m).node=mean(Nodes3D(Bricks(b,5:6),:),1);
                     elem(m).vertf(1:4,:)=Nodes3D(Bricks(b,5:8),:);
                     elem(m).vertf(5,:)=Nodes3D(Bricks(b,5),:);
                     elem(m).vertf(6,:)=Nodes3D(Bricks(b,6),:);
                     elem(m).vertf(7,:)=Nodes3D(Bricks(b,2),:);
                     elem(m).vertf(8,:)=Nodes3D(Bricks(b,1),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);

                     if i==1
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(Nt,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;

                     elseif i==Nt
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;

                     else
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;
                     end
                 else
                     b=bfun(i,j,h);
                     elem(m).type='cq';
                     elem(m).node=mean(Nodes3D(Bricks(b,5:8),:),1);
                     elem(m).vertf(1:4,:)=Nodes3D(Bricks(b,5:8),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);

                     if i==1
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(Nt,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;

                     elseif i==Nt
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;

                     else
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;
                     end
                 end
                 elem(m).Af=reshape(Aft(j,h,:),1,6);
                 elem(m).Ac=reshape(Act(j,h,:),1,6);
                 elem(m).V=V(j,h);
                 elem(m).dz_local=Geom_r(j,h);

                else
                 elem(m).face=[];
                 if j==Ntheta
                     elem(m).face=face_lateral(i);
                     elem(m).type='cq';
                     b=bfun(i,j,h);
                     elem(m).node=mean(Nodes3D(Bricks(b,[1,2,5,6]),:),1);
                     elem(m).vertf(1,:)=Nodes3D(Bricks(b,1),:);
                     elem(m).vertf(2,:)=Nodes3D(Bricks(b,2),:);
                     elem(m).vertf(3,:)=Nodes3D(Bricks(b,6),:);
                     elem(m).vertf(4,:)=Nodes3D(Bricks(b,5),:);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);

                     if i==1
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(Nt,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;

                     elseif i==Nt
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;

                     else
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;
                     end

                 else
                     b=bfun(i,j,h);
                     elem(m).type='i';
                     elem(m).node=mean(Nodes3D(Bricks(b,:),:),1);
                     elem(m).node_diff=mean(Nodes3D(Bricks(b,:)),1);

                     if i==1
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(Nt,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;

                     elseif i==Nt
                     Con(m,k(i-Nt+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;

                     else
                     Con(m,k(i+1,j,h))=1;
                     Con(m,k(i,j+1,h))=2;
                     Con(m,k(i,j,h+1))=3;
                     Con(m,k(i-1,j,h))=4;
                     if j>1
                     Con(m,k(i,j-1,h))=5;
                     end
                     Con(m,k(i,j,h-1))=6;
                     end

                 end
                 elem(m).Af=reshape(Afl(j,h,:),1,6);
                 elem(m).Ac=reshape(Acl(j,h,:),1,6);
                 elem(m).V=V(j,h);
                 elem(m).dz_local=Geom_r(j,h);
                end
            end

        end
    end
end

end
