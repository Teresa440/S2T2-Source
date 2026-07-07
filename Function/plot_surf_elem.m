function  plot_surf_elem(sat,plot_normals,n_fig)


surfaces=sat.geom.surfaces;
figure(n_fig)

for i=1:1:length(surfaces)
    for j=1:1:length(surfaces(i).elem)
        
        patch('x',surfaces(i).elem(j).vertf(:,1),'y',surfaces(i).elem(j).vertf(:,2),'z',surfaces(i).elem(j).vertf(:,3),'facecolor','y','edgecolor','r','FaceAlpha',0.1);
        hold on
    end
end
ax=linspace(0,400,200);
ax_z=zeros(size(ax));

plot3(ax,ax_z,ax_z,'r',LineWidth=3);
hold on
plot3(ax_z,ax,ax_z,'g',LineWidth=3);
hold on
plot3(ax_z,ax_z,ax,'b',LineWidth=3);
view(3)
daspect([1 1 1]);
xlabel('X'); ylabel('Y'); zlabel('Z');

% normals
if(plot_normals)
    for i = 1:1:length(sat.geom.surfaces)
        coord = [sat.geom.surfaces(i).center; sat.geom.surfaces(i).center + sat.geom.surfaces(i).norm*15];
        plot3(coord(:,1),coord(:,2),coord(:,3),'b',LineWidth=2);
        hold on
    end
end
end