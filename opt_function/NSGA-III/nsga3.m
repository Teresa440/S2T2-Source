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

function [pos_pareto, cost_pareto, pos_all, cost_all, pop] = nsga3(CostFunction, VarMin, VarMax, options)

if size(VarMin,2) ~= size(VarMax,2)
    error("VarMin and VarMax must have the same size")
end
% Convert to column vector
if size(VarMin,2) > size(VarMin,1)
    VarMin = VarMin';
    VarMax = VarMax';
    CostFunction = @(x) CostFunction(x');
end

VarSize = size(VarMin); % Size of Decision Variables Matrix

% Number of Objective Functions
test_val = CostFunction(rand(VarSize).*(VarMax - VarMin) + VarMin);
nObj = numel(test_val);
% Convert CostFunction output to be a column vector
if size(test_val,2) > size(test_val,1)
    CostFunction = @(x) CostFunction(x)';
end

%% NSGA-III Parameters

if ~exist('options','var')
    if nObj >= 2 && nObj <= 15
        suggested_ref_points = [0 100 12 8 6 4 4 3 3 3 2 2 2 2 2];
        options.nDivision = suggested_ref_points(N_objectives); % number of reference point in each axis of the nobjs dimensions
    end
    options.nPop = 100;  % Population Size
    options.pCrossover = 0.5;       % Crossover Percentage
    options.pMutation = 0.5;       % Mutation Percentage
    options.mu = 0.02;     % Mutation Rate
    options.MaxIt = 100;  % Maximum Number of Iterations
    options.plot_results = false;
    options.verbose = true;
end

nPop = options.nPop;

nCrossover = 2*round(options.pCrossover*options.nPop/2); % Number of Parnets (Offsprings)

nMutation = round(options.pMutation*options.nPop); % Number of Mutants

sigma = 0.1.*(VarMax-VarMin); % Mutation Step Size

% Generating Reference Points
Zr = GenerateReferencePoints(nObj, options.nDivision);

%% Collect Parameters

params.nPop = nPop;
params.Zr = Zr;
params.nZr = size(Zr,2);
params.zmin = [];
params.zmax = [];
params.smin = [];

%% Initialization

if options.verbose
    disp('Staring NSGA-III ...');
end

empty_individual.Position = [];
empty_individual.Cost = [];
empty_individual.Rank = [];
empty_individual.DominationSet = [];
empty_individual.DominatedCount = [];
empty_individual.NormalizedCost = [];
empty_individual.AssociatedRef = [];
empty_individual.DistanceToAssociatedRef = [];

pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    pop(i).Position = rand(VarSize).*(VarMax - VarMin) + VarMin;
    pop(i).Cost = CostFunction(pop(i).Position);
end

% Sort Population and Perform Selection
[pop, ~, params] = SortAndSelectPopulation(pop, params);


%% NSGA-III Main Loop

for it = 1:options.MaxIt
 
    % Crossover
    popc = repmat(empty_individual, nCrossover/2, 2);
    for k = 1:nCrossover/2

        i1 = randi([1 nPop]);
        p1 = pop(i1);

        i2 = randi([1 nPop]);
        p2 = pop(i2);

        [popc(k, 1).Position, popc(k, 2).Position] = Crossover_nsga3(p1.Position, p2.Position);

        popc(k, 1).Cost = CostFunction(popc(k, 1).Position);
        popc(k, 2).Cost = CostFunction(popc(k, 2).Position);

    end
    popc = popc(:);

    % Mutation
    popm = repmat(empty_individual, nMutation, 1);
    for k = 1:nMutation

        i = randi([1 nPop]);
        p = pop(i);

        popm(k).Position = Mutate(p.Position, options.mu, sigma, VarMin, VarMax);

        popm(k).Cost = CostFunction(popm(k).Position);

    end

    % Merge
    pop = [pop
           popc
           popm]; %#ok
    
    % Sort Population and Perform Selection
    [pop, F, params] = SortAndSelectPopulation(pop, params);
    
    % Store F1
    F1 = pop(F{1});

    % Show Iteration Information
    if options.verbose
        disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1))]);
    end

    % Plot F1 Costs
    if options.plot_results
        PlotCosts_ngsa3(F1,options.ax,it,options.MaxIt );
    end
 
end

%% Results

if options.verbose
    disp(['Final Iteration: Number of F1 Members = ' num2str(numel(F1))]);
    disp('Optimization Terminated.');
end

pos_pareto = [F1.Position]';
cost_pareto = [F1.Cost]';

pos_all = [pop.Position]';
cost_all = [pop.Cost]';

end


