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

function [pos_pareto, cost_pareto, pos_all, cost_all, pop] = moead(CostFunction, VarMin, VarMax, options)

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

%% MOEA/D Settings

if ~exist('options','var')
    options.nPop = 100; % Population Size (Number of Sub-Problems)
    options.nArchive=100;
    options.MaxIt=100;  % Maximum Number of Iterations
    options.gamma = 0.5;
    options.plot_results = false;
    options.verbose = true;
end

T=max(ceil(0.15*options.nPop),2);    % Number of Neighbors
T=min(max(T,2),15);

crossover_params.gamma = options.gamma;
crossover_params.VarMin=VarMin;
crossover_params.VarMax=VarMax;

% nFunEval = options.nPop.*options.MaxIt;

%% Initialization

% Create Sub-problems
sp=CreateSubProblems(nObj,options.nPop,T);

% Empty Individual
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.g=[];
empty_individual.IsDominated=[];

% Initialize Goal Point
%z=inf(nObj,1);
z=zeros(nObj,1);

% Create Initial Population
pop=repmat(empty_individual,options.nPop,1);
for i=1:options.nPop
    pop(i).Position=rand(VarSize).*(VarMax - VarMin) + VarMin;
    pop(i).Cost=CostFunction(pop(i).Position);
    z=min(z,pop(i).Cost);
end

for i=1:options.nPop
    pop(i).g=DecomposedCost(pop(i),z,sp(i).lambda);
end

% Determine Population Domination Status
pop=DetermineDomination_moead(pop);

% Initialize Estimated Pareto Front
EP=pop(~[pop.IsDominated]);

%% Main Loop

for it=1:options.MaxIt
    for i=1:options.nPop
        
        % Reproduction (Crossover)
        % K=randsample(T,2);
        K = randi(T,2,1);
        while K(1) == K(2)
            K = randi(T,2,1);
        end
        
        j1=sp(i).Neighbors(K(1));
        p1=pop(j1);
        
        j2=sp(i).Neighbors(K(2));
        p2=pop(j2);
        
        y=empty_individual;
        y.Position=Crossover_moead(p1.Position,p2.Position,crossover_params);
        
        y.Cost=CostFunction(y.Position);
        
        z=min(z,y.Cost);
        
        for j=sp(i).Neighbors
            y.g=DecomposedCost(y,z,sp(j).lambda);
            if y.g<=pop(j).g
                pop(j)=y;
            end
        end
        
    end
    
    % Determine Population Domination Status
	pop=DetermineDomination_moead(pop);
    
    ndpop=pop(~[pop.IsDominated]);
    
    EP=[EP
        ndpop]; %#ok
    
    EP=DetermineDomination_moead(EP);
    EP=EP(~[EP.IsDominated]);
    
    % if numel(EP)>nArchive
    %     Extra=numel(EP)-nArchive;
    %     ToBeDeleted=randsample(numel(EP),Extra);
    %     EP(ToBeDeleted)=[];
    % end
    while numel(EP) > options.nArchive
        ToBeDeleted=randi(numel(EP));
        EP(ToBeDeleted)=[];
    end

    if options.plot_results
        % Plot EP
        % figure(1);
        PlotCosts_moead(EP,options.ax,it,options.MaxIt);
        pause(0.01);
    end
    
    if options.verbose
        % Display Iteration Information
        disp(['Iteration ' num2str(it) ': Number of Pareto Solutions = ' num2str(numel(EP))]);
    end
    
end

%% Reults

EPC=[EP.Cost];
if options.verbose
    disp(' ');
    for j=1:nObj
        
        disp(['Objective #' num2str(j) ':']);
        disp(['      Min = ' num2str(min(EPC(j,:)))]);
        disp(['      Max = ' num2str(max(EPC(j,:)))]);
        disp(['    Range = ' num2str(max(EPC(j,:))-min(EPC(j,:)))]);
        disp(['    St.D. = ' num2str(std(EPC(j,:)))]);
        disp(['     Mean = ' num2str(mean(EPC(j,:)))]);
        disp(' ');
        
    end
end

pos_all = [EP.Position]';
cost_all = [EP.Cost]';
is_dominated = [EP.IsDominated]';

pos_pareto = pos_all(is_dominated == 0,:);
cost_pareto = cost_all(is_dominated == 0,:);

end


