% 
% Copyright (c) 2016, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
% 
% Project Code: YPEA126
% Project Title: Non-dominated Sorting Genetic Algorithm III (NSGA-III)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Implemented by: S. Mostapha Kalami Heris, PhD (member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
% 
% Base Reference Paper:
% K. Deb and H. Jain, "An Evolutionary Many-Objective Optimization Algorithm 
% Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving
% Problems With Box Constraints,"
% in IEEE Transactions on Evolutionary Computation,
% vol. 18, no. 4, pp. 577-601, Aug. 2014.
% 
% Reference Papaer URL: http://doi.org/10.1109/TEVC.2013.2281535
% 

function PlotCosts_ngsa3(pop,graph,it,MaxIt)

    Costs=[pop.Cost];
    
    if size(Costs,1) == 2
    plot(graph,Costs(1,:),Costs(2,:),'ok');
    elseif size(Costs,1) > 2
        plot3(graph,Costs(1,:),Costs(2,:),Costs(3,:),'ok');
        zlabel(graph,'3^{rd} Objective');
    end
    xlabel(graph,'1^{st} Objective');
    ylabel(graph,'2^{nd} Objective');
    title(graph,"Iteration " + string(it) + "/" + string(MaxIt))
    grid(graph,"on")
    hold(graph,"off")
    drawnow

end