%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA124
% Project Title: Implementation of MOEA/D
% Muti-Objective Evolutionary Algorithm based on Decomposition
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function PlotCosts_moead(EP,graph,it,MaxIt)

    EPC=[EP.Cost];
    dominated = [EP.IsDominated];
    EPC_dominated = EPC(:,dominated == 1);
    EPC_pareto = EPC(:,dominated == 0);
    if size(EPC,1) == 2
    plot(graph,EPC_pareto(1,:),EPC_pareto(2,:),'ok');
    hold(graph,"on")
    plot(graph,EPC_dominated(1,:),EPC_dominated(2,:),'or');
    elseif size(EPC,1) > 2
        plot3(graph,EPC_pareto(1,:),EPC_pareto(2,:),EPC_pareto(3,:),'ok');
        hold(graph,"on")
        plot3(graph,EPC_dominated(1,:),EPC_dominated(2,:),EPC_dominated(3,:),'or');
        zlabel(graph,'3^{rd} Objective');
    end
    xlabel(graph,'1^{st} Objective');
    ylabel(graph,'2^{nd} Objective');
    title(graph,"Generation " + string(it) + "/" + string(MaxIt))
    grid(graph,"on")
    hold(graph,"off")
    drawnow
end