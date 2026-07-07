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

% function y=Mutate(x,mu,sigma)
% 
%     nVar=numel(x);
%     
%     nMu=ceil(mu*nVar);
% 
%     j=randsample(nVar,nMu);
% 
%     y=x;
%     
%     y(j)=x(j)+sigma*randn(size(j));
% 
% end

function y=Mutate(x,mu,sigma,VarMin,VarMax)

    nVar=numel(x);
    
    j = rand(nVar,1) < mu;

    if sum(j) == 0
        j(randi(nVar)) = 1;
    end

    y=x;
    
    y(j)=x(j)+sigma(j).*(randn(sum(j),1));

    y = regularize_mutant(y, VarMin, VarMax);

    function new_pt = regularize_mutant(new_pt, lb, ub)
        %REGULARIZE_MUTANT: regularize infeasible mutant, honoring lower and upper bounds
    
        cont = 0;
        under_bound = new_pt < lb;
        while sum(new_pt < lb) ~= 0
            cont = cont + 1;
            new_pt(under_bound) = new_pt(under_bound) + ub(under_bound) - lb(under_bound);
            if cont == 2
                new_pt(under_bound) = ub(under_bound);
                break
            end
        end
        cont = 0;
        over_bound = new_pt > ub;
        while sum(new_pt > ub) ~= 0
            cont = cont + 1;
            new_pt(over_bound) = new_pt(over_bound) - ub(over_bound) + lb(over_bound);
            if cont == 2
                new_pt(over_bound) = lb(over_bound);
            break
            end
        end
    end

end