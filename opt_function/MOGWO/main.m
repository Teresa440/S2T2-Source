clear all
clc

TestProblem='UF1';
nVar=10;

fobj = cec09(TestProblem);

xrange = xboundary(TestProblem, nVar);

% Lower bound and upper bound
lb=xrange(:,1);
ub=xrange(:,2);

% options structure
options = struct;
options.GreyWolves_num = 100;
options.MaxIt = 100; % Maximum Number of Iterations
options.Archive_size = 100; % Repository Size
options.alpha=0.1; % Grid Inflation Parameter
options.nGrid=10; % Number of Grids per each Dimension
options.beta=4; % Leader Selection Pressure Parameter
options.gamma=2; % Extra (to be deleted) Repository Member Selection Pressure
options.plot_results = true;
options.verbose = true;

[pos_pareto, cost_pareto, pos_all, cost_all] = MOGWO(fobj, lb, ub, options);

figure
plot(cost_pareto(:,1),cost_pareto(:,2),'xk');
