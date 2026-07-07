function state = custom_gamultiobj_plot(options,state,~,graph)

n_obj = size(state.Score,2); % number of objectives
if n_obj == 2
    plot(graph,state.Score(:,1),state.Score(:,2),'ok')
%     hold(graph,"on")
%     plot(graph,state.Best(:,1),state.Best(:,2),'ok')
elseif n_obj > 2
    plot3(graph,state.Score(:,1),state.Score(:,2),state.Score(:,3),'ok') 
%     hold(graph,"on")
%     plot3(graph,state.Best(:,1),state.Best(:,2),state.Best(:,3),'ok') 
    zlabel(graph,'3^{rd} Objective'); 
end
xlabel(graph,'1^{st} Objective'); ylabel(graph,'2^{nd} Objective');
title(graph,"Generation " + string(state.Generation) + "/" + string(options.Generations))
hold(graph,"off")
grid(graph,"on")

end