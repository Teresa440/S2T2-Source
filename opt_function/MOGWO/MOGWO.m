%___________________________________________________________________%
%  Multi-Objective Grey Wolf Optimizer (MOGWO)                      %
%  Source codes demo version 1.0                                    %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%    S. Mirjalili, S. Saremi, S. M. Mirjalili, L. Coelho,           %
%    Multi-objective grey wolf optimizer: A novel algorithm for     %
%    multi-criterion optimization, Expert Systems with Applications,%
%    in press, DOI: http://dx.doi.org/10.1016/j.eswa.2015.10.039    %
%                                                                   %
%___________________________________________________________________%

% I acknowledge that this version of MOGWO has been written using
% a large portion of the following code:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MATLAB Code for                                                  %
%                                                                   %
%  Multi-Objective Particle Swarm Optimization (MOPSO)              %
%  Version 1.0 - Feb. 2011                                          %
%                                                                   %
%  According to:                                                    %
%  Carlos A. Coello Coello et al.,                                  %
%  "Handling Multiple Objectives with Particle Swarm Optimization," %
%  IEEE Transactions on Evolutionary Computation, Vol. 8, No. 3,    %
%  pp. 256-279, June 2004.                                          %
%                                                                   %
%  Developed Using MATLAB R2009b (Version 7.9)                      %
%                                                                   %
%  Programmed By: S. Mostapha Kalami Heris                          %
%                                                                   %
%         e-Mail: sm.kalami@gmail.com                               %
%                 kalami@ee.kntu.ac.ir                              %
%                                                                   %
%       Homepage: http://www.kalami.ir                              %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos_pareto, cost_pareto, pos_all, cost_all] = MOGWO(CostFunction, lb, ub, options)

% Convert input to row vector
if size(lb,2) < size(lb,1)
    lb = lb';
    ub = ub';
    CostFunction = @(x) CostFunction(x');
end

VarSize = size(lb); % Size of Decision Variables Matrix
nVar = numel(lb);
% Number of Objective Functions
test_val = CostFunction(rand(VarSize).*(ub - lb) + lb);
% Convert CostFunction output to be a row vector
if size(test_val,2) < size(test_val,1)
    CostFunction = @(x) CostFunction(x)';
end

%basic settings
if ~exist('options','var')
    options = struct;
    options.GreyWolves_num = 100;
    options.MaxIt = 100; % Maximum Number of Iterations
    options.Archive_size = 100; % Repository Size
    options.alpha=0.1; % Grid Inflation Parameter
    options.nGrid=10; % Number of Grids per each Dimension
    options.beta=4; % Leader Selection Pressure Parameter
    options.gamma=2; % Extra (to be deleted) Repository Member Selection Pressure
    options.plot_results = false;
    options.verbose = true;
end
GreyWolves_num = options.GreyWolves_num;
MaxIt = options.MaxIt;
Archive_size = options.Archive_size;
alpha = options.alpha;
nGrid = options.nGrid;
beta = options.beta;
gamma = options.gamma;
% f_eval = options.GreyWolves_num*options.MaxIt;

% Initialization

GreyWolves=CreateEmptyParticle(GreyWolves_num);


for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    GreyWolves(i).Position=zeros(1,nVar);
    for j=1:nVar
        % GreyWolves(i).Position(1,j)=unifrnd(lb(j),ub(j),1);
        GreyWolves(i).Position(1,j)=rand(1).*(ub(j) - lb(j)) + lb(j);
    end
    % GreyWolves(i).Cost=CostFunction(GreyWolves(i).Position')';
    GreyWolves(i).Cost=CostFunction(GreyWolves(i).Position);
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end

GreyWolves=DetermineDomination_MOGWO(GreyWolves);

Archive=GetNonDominatedParticles(GreyWolves);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

% MOGWO main loop

for it=1:MaxIt
    a=2-it*((2)/MaxIt);
    for i=1:GreyWolves_num
        
        clear rep2
        clear rep3
        
        % Choose the alpha, beta, and delta grey wolves
        Delta=SelectLeader(Archive,beta);
        Beta=SelectLeader(Archive,beta);
        Alpha=SelectLeader(Archive,beta);
        
        % If there are less than three solutions in the least crowded
        % hypercube, the second least crowded hypercube is also found
        % to choose other leaders from.
        if size(Archive,1)>1
            counter=0;
            for newi=1:size(Archive,1)
                if sum(Delta.Position~=Archive(newi).Position)~=0
                    counter=counter+1;
                    rep2(counter,1)=Archive(newi);
                end
            end
            Beta=SelectLeader(rep2,beta);
        end
        
        % This scenario is the same if the second least crowded hypercube
        % has one solution, so the delta leader should be chosen from the
        % third least crowded hypercube.
        if size(Archive,1)>2
            counter=0;
            for newi=1:size(rep2,1)
                if sum(Beta.Position~=rep2(newi).Position)~=0
                    counter=counter+1;
                    rep3(counter,1)=rep2(newi);
                end
            end
            Alpha=SelectLeader(rep3,beta);
        end
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Delta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand(1, nVar)-a;
        % Eq.(3.8) in the paper
        X1=Delta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Beta.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.9) in the paper
        X2=Beta.Position-A.*abs(D);
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nVar);
        % Eq.(3.1) in the paper
        D=abs(c.*Alpha.Position-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand()-a;
        % Eq.(3.10) in the paper
        X3=Alpha.Position-A.*abs(D);
        
        % Eq.(3.11) in the paper
        GreyWolves(i).Position=(X1+X2+X3)./3;
        
        % % Boundary checking
        % GreyWolves(i).Position=min(max(GreyWolves(i).Position,lb),ub);
        %
        % if ~all(GreyWolves(i).Position <= ub & GreyWolves(i).Position >= lb)
        %     error("Error: out of bounds")
        % end
        GreyWolves(i).Position = regularize_mutant(GreyWolves(i).Position, lb, ub);

        % GreyWolves(i).Cost=CostFunction(GreyWolves(i).Position')';
        GreyWolves(i).Cost=CostFunction(GreyWolves(i).Position);
    end
    
    GreyWolves=DetermineDomination_MOGWO(GreyWolves);
    non_dominated_wolves=GetNonDominatedParticles(GreyWolves);
    
    Archive=[Archive
        non_dominated_wolves];
    
    Archive=DetermineDomination_MOGWO(Archive);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alpha);
        
    end
    
    if options.verbose==1
        disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    end
    % save results
    
    % Results
    
    costs=GetCosts(GreyWolves);
    Archive_costs=GetCosts(Archive);
    
    if options.plot_results==1
        hold(options.ax,"off")
        if size(costs,1) == 2
            plot(options.ax,costs(1,:),costs(2,:),'or');        
            hold(options.ax,"on")
            plot(options.ax,Archive_costs(1,:),Archive_costs(2,:),'ok');
            % legend(options.ax,'Grey wolves','Non-dominated solutions');        
        elseif size(costs,1) > 2
            plot3(options.ax,costs(1,:),costs(2,:),costs(3,:),'or');        
            hold(options.ax,"on")
            plot3(options.ax,Archive_costs(1,:),Archive_costs(2,:),Archive_costs(3,:),'ok');
            zlabel(options.ax,'3^{rd} Objective'); 
        end
        xlabel(options.ax,'1^{st} Objective'); ylabel(options.ax,'2^{nd} Objective');
        title(options.ax,"Iteration " + string(it) + "/" + string(MaxIt))
        grid(options.ax,"on")
        drawnow
    end
    
end

pos_all = zeros(length(GreyWolves), length(GreyWolves(1).Position));
cost_all = zeros(length(GreyWolves), length(GreyWolves(1).Cost));
for i = 1:length(GreyWolves)
    pos_all(i,:) = GreyWolves(i).Position;
    cost_all(i,:) = GreyWolves(i).Cost;
end

pos_pareto = zeros(length(Archive), length(Archive(1).Position));
cost_pareto = zeros(length(Archive), length(Archive(1).Cost));
for i = 1:length(Archive)
    pos_pareto(i,:) = Archive(i).Position;
    cost_pareto(i,:) = Archive(i).Cost;
end

end

%% Local Functions

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

