
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  The source code of the reference vector guided evolutionary algorithm (RVEA)
%
%  See the details of RVEA in the following paper:
%
%  R. Cheng, Y. Jin, M. Olhofer and B. Sendhoff, 
%  A Reference Vector Guided Evolutionary Algorithm for Many-objective Optimization,
%  IEEE Transactions on Evolutionary Computation, 2016
%
%  The source code RVEA is implemented by Ran Cheng 
%
%  If you have any questions about the code, please contact: 
%  
%  Ran Cheng at ranchengcn@gmail.com
%  Prof. Yaochu Jin at yaochu.jin@surrey.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos_pareto, cost_pareto, pos_all, cost_all] = RVEA(CostFunction, VarMin, VarMax, options)

nEval = 0;

% Convert input to row vector (need to be fixed)
if size(VarMin,2) < size(VarMin,1)
    VarMin = VarMin';
    VarMax = VarMax';
    CostFunction = @(x) CostFunction(x'); 
end

VarSize = size(VarMin); % Size of Decision Variables Matrix
% Number of Objective Functions
test_val = CostFunction(rand(VarSize).*(VarMax - VarMin) + VarMin);
M = numel(test_val);
% Convert CostFunction output to be a row vector
if size(test_val,2) < size(test_val,1)
    CostFunction = @(x) CostFunction(x)';
    nEval = nEval + 1;
end


%basic settings
if ~exist('options','var')
    options = struct;
    p1_v = [99 13  7  5  4  3  3  2  3];
    p2_v = [ 0  0  0  0  1  2  2  2  2];
    options.p1 = p1_v(M-1);
    options.p2 = p2_v(M-1);
    options.alpha = 2.0;
    options.FE = 0;
    options.MaxFE = 20000;
    options.ref_vect_adapt_step = 10;
    options.plot_results = true;
    options.verbose = true;
end
p1 = options.p1;
p2 = options.p2;
alpha = options.alpha;
FE = options.FE;
% Evaluations = Generations*N;
% N = nchoosek(p1+M-1,M-1);
% if p2 > 0
% N = N + nchoosek(p2+M-1,M-1);
% end

%reference vector initialization
[N,Vs] = F_weight(p1,p2,M);
Vs(Vs==0) = 0.000000001;
for i = 1:N
    Vs(i,:) = Vs(i,:)./norm(Vs(i,:));
end
V = Vs;

%calculat neighboring angle for angle normalization
cosineVV = V*V';
[scosineVV, ~] = sort(cosineVV, 2, 'descend');
acosVV = acos(scosineVV(:,2));
refV = (acosVV);

%population initialization
Population = rand(N,length(VarMin));
Population = Population.*repmat(VarMax,N,1)+(1-Population).*repmat(VarMin,N,1);
Boundary = [VarMax;VarMin];
Coding   = 'Real';

FunctionValue = zeros(N,M);
for vc = 1:size(Population,1)
    FunctionValue(vc,:) = CostFunction(Population(vc,:));
    nEval = nEval + 1;
end
% FunctionValue = P_objective('value',Problem,M,Population);

Gene = 0;
while FE < options.MaxFE
    %random mating and reproduction
    [MatingPool] = F_mating(Population);
    Offspring = P_generator(MatingPool,Boundary,Coding,N);  
    FE = FE + size(Offspring, 1);
    Population = [Population; Offspring];
    FunctionValueAppend = zeros(size(Offspring,1),M);
    for vc = 1:size(Offspring,1)
        FunctionValueAppend(vc,:) = CostFunction(Offspring(vc,:));
        nEval = nEval + 1;
    end
    FunctionValue = [FunctionValue; FunctionValueAppend];
    
    %APD based selection
    theta0 =  (FE/options.MaxFE)^alpha*(M);
    [Selection] = F_select(FunctionValue,V, theta0, refV);
    Population = Population(Selection,:);
    FunctionValue = FunctionValue(Selection,:);

    %reference vector adaption
    if(mod(Gene, options.ref_vect_adapt_step) == 0)
        %update the reference vectors
        Zmin = min(FunctionValue,[],1);	
        Zmax = max(FunctionValue,[],1);	
        V = Vs;
        V = V.*repmat((Zmax - Zmin)*1.0,N,1);
        for i = 1:N
            V(i,:) = V(i,:)./norm(V(i,:));
        end
        %update the neighborning angle value for angle normalization
        cosineVV = V*V';
        [scosineVV, ~] = sort(cosineVV, 2, 'descend');
        acosVV = acos(scosineVV(:,2));
        refV = (acosVV); 
    end

    if options.verbose
        fprintf('Generation: %d\tPopulation size = %d\tFunction Evaluation = %d\n',Gene,size(Population,1),FE);
    end

    if options.plot_results
        NonDominated = P_sort(FunctionValue,'first')==1;
        if size(FunctionValue,2) == 2
            plot(options.ax,FunctionValue(NonDominated,1),FunctionValue(NonDominated,2),'ok')
            hold(options.ax,"on")
            plot(options.ax,FunctionValue(~NonDominated,1),FunctionValue(~NonDominated,2),'or')
        elseif size(FunctionValue,2) > 2
            plot3(options.ax,FunctionValue(NonDominated,1),FunctionValue(NonDominated,2),FunctionValue(NonDominated,3),'ok')
            hold(options.ax,"on")
            plot3(options.ax,FunctionValue(~NonDominated,1),FunctionValue(~NonDominated,2),FunctionValue(~NonDominated,3),'or')
            zlabel(options.ax,'3^{rd} Objective');
        end
        xlabel(options.ax,'1^{st} Objective');
        ylabel(options.ax,'2^{nd} Objective');
        title(options.ax,"Function evaluation " + string(FE) + "/" + string(options.MaxFE))
        grid(options.ax,"on")
        hold(options.ax,"off")
        drawnow
    end
    Gene = Gene + 1;
end

% FunctionValue = P_objective('value',Problem,M,Population);
% for vc = 1:size(Population,1)
%     FunctionValue(vc,:) = CostFunction(Population(vc,:));
%     nEval = nEval + 1;
% end
NonDominated = P_sort(FunctionValue,'first')==1;
PopulationFront = Population(NonDominated,:);
FunctionValueFront = FunctionValue(NonDominated,:);

pos_pareto = PopulationFront;
cost_pareto = FunctionValueFront;
pos_all = Population;
cost_all = FunctionValue;

if options.verbose
    disp("Total Function Evaluation: " + string(nEval))
end

end


