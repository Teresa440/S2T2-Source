function plot_GMM4_presentation(sat,graph,vis)

if exist("vis","var")
    alpha_structure = 1-vis;
else
    alpha_structure = 0.1;
end

lw = 0.2;

for i=1:1:length(sat.geom.ext.face)
    surf(graph,sat.geom.ext.face(i).gridX,sat.geom.ext.face(i).gridY,sat.geom.ext.face(i).gridZ,'EdgeColor',[0, 38, 79]./256,'FaceColor','c',...
        'FaceAlpha',alpha_structure,'LineStyle','-','LineWidth',lw);
    hold(graph,'on')
end

for i=1:1:sat.geom.Nsp
    cc = sat.geom.sp(i).color;
    sp=surf(graph,sat.geom.sp(i).face.gridX,sat.geom.sp(i).face.gridY,sat.geom.sp(i).face.gridZ);
    sp.FaceColor=string(cc);
    sp.EdgeColor=[0, 38, 79]./256;
    sp.LineStyle='-';
    sp.LineWidth=lw;
    sp.FaceAlpha=0.7;
    hold(graph,'on')
end

for i=1:1:sat.geom.Nb
    cc = sat.geom.board(i).color;
    b=surf(graph,sat.geom.board(i).face.gridX,sat.geom.board(i).face.gridY,sat.geom.board(i).face.gridZ);
    b.FaceColor=string(cc);
    b.EdgeColor=[0, 38, 79]./256;
    b.LineStyle='-';
    b.LineWidth=lw;
    b.FaceAlpha=0.7;
    hold(graph,'on')
end

for j=1:1:sat.geom.NP
    cc = sat.geom.parall(j).color;
    for i=1:6        
        p=surf(graph,sat.geom.parall(j).face(i).gridX,...
             sat.geom.parall(j).face(i).gridY,...
             sat.geom.parall(j).face(i).gridZ);
        p.FaceColor=string(cc);
        p.EdgeColor=[0, 38, 79]./256;
        p.LineStyle='-';
        p.LineWidth=lw;
        p.FaceAlpha=0.7;
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
    cc = sat.geom.cyl(j).color;
    for i=1:1:length(elem)
        
        plot3(graph,elem(i).node(1),elem(i).node(2),elem(i).node(3),'-o','Color','b','MarkerSize',4,...
            'MarkerFaceColor','w')
        hold(graph,'on')
        if isempty(elem(i).vertf)==0
            if strcmp(elem(i).type,'s')==1
              p(i*2)=patch(graph,elem(i).vertf(1:4,1),elem(i).vertf(1:4,2),elem(i).vertf(1:4,3),'y');
              p(i*2).FaceColor = string(cc);
              p(i*2).EdgeColor = [0, 38, 79]./256;
              p(i*2).LineStyle='-';
              p(i*2).LineStyle=lw;
              p(i*2).FaceAlpha = 0.7;
              hold(graph,'on')
              p(i*2+1)=patch(graph,elem(i).vertf(5:8,1),elem(i).vertf(5:8,2),elem(i).vertf(5:8,3),'y');
              p(i*2+1).FaceColor = string(cc);
              p(i*2+1).EdgeColor = [0, 38, 79]./256;
              p(i*2+1).LineStyle='-';
              p(i*2+1).LineStyle=lw;
              p(i*2+1).FaceAlpha = 0.7;
              hold(graph,'on')
            else
                p(i*2)=patch(graph,elem(i).vertf(:,1),elem(i).vertf(:,2),elem(i).vertf(:,3),'y');
                p(i*2).FaceColor = string(cc);
                p(i*2).EdgeColor = [0, 38, 79]./256;
                p(i*2).LineStyle='-';
                p(i*2).LineStyle=lw;
                p(i*2).FaceAlpha = 0.7;
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

set(gca, 'Color','w', 'XColor',[0, 38, 79]./256, 'YColor',[0, 38, 79]./256, 'ZColor',[0, 38, 79]./256)

end
