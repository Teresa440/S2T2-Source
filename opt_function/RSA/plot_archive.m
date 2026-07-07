function plot_archive(values, dominance, ref_points, ax,  cont, T, Tmax, Tmin)
%PLOT_ARCHIVE: plots in the objective space the points of the archive

if size(values,2) == 2
    plot(ax,values(dominance==2,1),values(dominance==2,2),"ok")
    hold(ax,"on")
    plot(ax,values(dominance==1,1),values(dominance==1,2),"or")
    % for k = 1:size(ref_points,1)
        % plot([0 ref_points(k,1)*10], [0 ref_points(k,2)*10], ":k","MarkerSize",0.2,"LineWidth",0.2)
    % end
    if isempty(values(dominance==1,1))
        % legend("Non-dominated points","Reference Lines")
        % legend("Non-dominated points")
    else
        % legend("Non-dominated points","Dominated points","Reference Lines")
        % legend("Non-dominated points","Dominated points")
    end
    xlabel(ax,'1^{st} Objective'); ylabel(ax,'2^{nd} Objective'); 
elseif size(values,2) > 2
    plot3(ax,values(dominance==2,1),values(dominance==2,2),values(dominance==2,3),"ok")
    hold(ax,"on")
    plot3(ax,values(dominance==1,1),values(dominance==1,2),values(dominance==1,3),"or")
    % for k = 1:size(ref_points,1)
       % plot3([0 ref_points(k,1)*10], [0 ref_points(k,2)*10], [0 ref_points(k,3)*10], ":k","MarkerSize",0.2,"LineWidth",0.2)
    % end
    if isempty(values(dominance==1,1))
        % legend("Non-dominated points","Reference Lines")
        % legend("Non-dominated points")
    else
        % legend("Non-dominated points","Dominated points","Reference Lines")
        % legend("Non-dominated points","Dominated points")
    end
    xlabel(ax,'1^{st} Objective'); ylabel(ax,'2^{nd} Objective'); zlabel(ax,'3^{rd} Objective'); 
else
    disp("Error: only 2 or 3 dimensions objective spaces can be plotted")
    return
end
title(ax,"Iteration " + string(cont) + "    Temperature " + string(Tmax) + " -> " + string(T) + " -> " + string(Tmin))
grid(ax,"on")
hold(ax,"off")
drawnow

end

