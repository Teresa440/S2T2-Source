function plot_GMM4_presentation(sat,graph,vis)

if exist("vis","var")
    alpha_structure = 1-vis;
else
    alpha_structure = 0;
end

lw = 0.2;
alpha2 = 1;
face_c = [0, 38, 79]./256;
edge_c = [229,240,255]./255;

for i=1:1:length(sat.geom.ext.face)
    surf(graph,sat.geom.ext.face(i).gridX,sat.geom.ext.face(i).gridY,sat.geom.ext.face(i).gridZ,'EdgeColor',edge_c,'FaceColor',face_c,...
        'FaceAlpha',alpha_structure,'LineStyle','-','LineWidth',lw);
    hold(graph,'on')
end

for i=1:1:sat.geom.Nsp
    sp=surf(graph,sat.geom.sp(i).face.gridX,sat.geom.sp(i).face.gridY,sat.geom.sp(i).face.gridZ);
    sp.FaceColor=face_c;
    sp.EdgeColor=edge_c;
    sp.LineStyle='-';
    sp.LineWidth=lw;
    sp.FaceAlpha=alpha2;
    hold(graph,'on')
end

for i=1:1:sat.geom.Nb
    b=surf(graph,sat.geom.board(i).face.gridX,sat.geom.board(i).face.gridY,sat.geom.board(i).face.gridZ);
    b.FaceColor=face_c
    b.EdgeColor=edge_c;
    b.LineStyle='-';
    b.LineWidth=lw;
    b.FaceAlpha=alpha2;
    hold(graph,'on')
end

for j=1:1:sat.geom.NP
    for i=1:6        
        p=surf(graph,sat.geom.parall(j).face(i).gridX,...
             sat.geom.parall(j).face(i).gridY,...
             sat.geom.parall(j).face(i).gridZ);
        p.FaceColor=face_c
        p.EdgeColor=edge_c;
        p.LineStyle='-';
        p.LineWidth=lw;
        p.FaceAlpha=alpha2;
         hold(graph,'on')       
         
    end
end

% for i=1:1:sat.geom.Nc
%     for j=1:1:length(sat.geom.cyl(i).face)
%         plot3(sat.geom.cyl(i).face(j).mesh(:,1),...
%               sat.geom.cyl(i).face(j).mesh(:,2),...
%               sat.geom.cyl(i).face(j).mesh(:,3));
%     end
% end

for j=1:1:sat.geom.Nc

    elem=sat.node.cyl(j).elements;
    for i=1:1:length(elem)
        
        plot3(graph,elem(i).node(1),elem(i).node(2),elem(i).node(3),'-o','Color','b','MarkerSize',4,...
            'MarkerFaceColor','w')
        hold(graph,'on')
        if isempty(elem(i).vertf)==0
            if strcmp(elem(i).type,'s')==1
              p(i*2)=patch(graph,elem(i).vertf(1:4,1),elem(i).vertf(1:4,2),elem(i).vertf(1:4,3),'y');
              p(i*2).FaceColor = face_c;
              p(i*2).EdgeColor = edge_c;
              p(i*2).LineStyle='-';
              p(i*2).LineStyle=lw;
              p(i*2).FaceAlpha = alpha2;
              hold(graph,'on')
              p(i*2+1)=patch(graph,elem(i).vertf(5:8,1),elem(i).vertf(5:8,2),elem(i).vertf(5:8,3),'y');
              p(i*2+1).FaceColor = face_c;
              p(i*2+1).EdgeColor = edge_c;
              p(i*2+1).LineStyle='-';
              p(i*2+1).LineStyle=lw;
              p(i*2+1).FaceAlpha = alpha2;
              hold(graph,'on')
            else
                p(i*2)=patch(graph,elem(i).vertf(:,1),elem(i).vertf(:,2),elem(i).vertf(:,3),'y');
                p(i*2).FaceColor = face_c;
                p(i*2).EdgeColor = edge_c;
                p(i*2).LineStyle='-';
                p(i*2).LineStyle=lw;
                p(i*2).FaceAlpha = alpha2;
                hold(graph,'on')
            end
        end
    end
end

% ax=linspace(0,400,200);
ax=[-20,sat.geom.ext.size(1)+20];
ay=[-20,sat.geom.ext.size(2)+20];
az=[-20,sat.geom.ext.size(3)+20];

ax_z=zeros(size(ax));

axis([ax, ay, az])

% plot3(graph,ax,ax_z,ax_z,'r',LineWidth=3);
% hold(graph,'on')
% plot3(graph,ax_z,ay,ax_z,'g',LineWidth=3);
% hold(graph,'on')
% plot3(graph,ax_z,ax_z,az,'b',LineWidth=3);

daspect(graph,[1 1 1]);
xlabel(graph,'X'); ylabel(graph,'Y'); zlabel(graph,'Z');

set(gca, 'Color',[0, 38, 79]./256, 'XColor',[229,240,255]./256, 'YColor',[229,240,255]./256, 'ZColor',[229,240,255]./256)

% axis off
grid off

drawnow

end
