





 plot_surf_elem(sat)



center=surf_for_MCRT(10).elem(1).center;

for i=5000:1:10000
    line=[center;VV(108).pend(i,:)];
    plot3(line(:,1),line(:,2),line(:,3))
    hold on 
end

% for i=1:1:length(surf_for_MCRT)
%     p1=surf_for_MCRT(i).center;
%     norm=surf_for_MCRT(i).norm;
%     p2=[p1(1)+norm(1)*30,p1(2)+norm(2)*30,p1(3)+norm(3)*30];
%     line=[p1;p2];
%     plot3(line(:,1),line(:,2),line(:,3),'g',LineWidth=3);
%     hold on
% end
% i=10;
% 
% for i=9:1:16
%      patch('x',surf_for_MCRT(i).vert(:,1),'y',surf_for_MCRT(i).vert(:,2),'z',surf_for_MCRT(i).vert(:,3),'facecolor','none','edgecolor','g');
%         hold on
% end
