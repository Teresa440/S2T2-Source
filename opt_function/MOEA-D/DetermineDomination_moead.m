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

function pop=DetermineDomination_moead(pop)

    nPop=numel(pop);

    for i=1:nPop
        pop(i).IsDominated=false;
    end
    
    for i=1:nPop
        for j=i+1:nPop
            if Dominates_moead(pop(i),pop(j))
                pop(j).IsDominated=true;
                
            elseif Dominates_moead(pop(j),pop(i))
                pop(i).IsDominated=true;
                
            end
        end
    end

end