function new_pt = perturb(archive, options, nvars, curr_point, curr_iter, lb, ub)
%PERTURB: perturbes the current solution to generate a new one
%   uses 4 different types of mutation, switching from one to the others

if size(curr_point,2) ~= nvars
    disp("Wrong size of curr_point");
end

archive_vars = archive.vars(archive.slot_status ~= 0,:);

if curr_iter/options.iter > options.switch_factor
    new_pt = differential_mut(archive_vars, nvars, curr_point, lb, ub);
    if options.flag2 == 1
        new_pt = polynomial_mut(new_pt, nvars, lb, ub);
    end
    if rand(1) < options.prob2
        new_pt = laplacian_mut(new_pt, nvars, lb, ub);
    end
else
    new_pt = SBX_mut(archive_vars, nvars, curr_point, lb, ub);
    if options.flag2 == 1
        new_pt = polynomial_mut(new_pt, nvars, lb, ub);
    end
    if rand(1) < options.prob2
        new_pt = laplacian_mut(new_pt, nvars, lb, ub);
    end
end

if sum(new_pt < lb) ~= 0 || sum(new_pt > ub) ~= 0 
    disp("Error: point out of bounds");
end

end

%% Local functions

function new_pt = SBX_mut(archive_vars, nvars, curr_point, lb, ub)
%SBX_mut: performs simulated binary crossover mutation

parent_2_id = randi(size(archive_vars, 1));
parent_2 = archive_vars(parent_2_id,:); % second parent

eta_c = 30; % tunable parameter (30 is suggested)
u = rand(1,nvars);
beta = (2*u).^(1/(eta_c + 1)).*(u <= 0.5) + (2*(1 - u)).^(-1/(eta_c + 1)).*(u > 0.5);

new_pt = 0.5*((1 + beta).*curr_point + (1 - beta).*parent_2);

new_pt = regularize_mutant(new_pt, lb, ub); % honor bounds

end


function new_pt = differential_mut(archive_vars, nvars, curr_point, lb, ub)
%DIFFERENTIAL_MUT: select a mutation base and 2 donors to generate a mutant
%v, then binomial crossover is performed between curr_point and v

b_id = randi(size(archive_vars, 1));
b = archive_vars(b_id,:);
x1_id = randi(size(archive_vars, 1));
x1 = archive_vars(x1_id,:); % first donor
x2_id = randi(size(archive_vars, 1));
x2 = archive_vars(x2_id,:); % second donor

F = 0.5; % tunable parameter (0.5 is suggested), could also be computed using Laplace distribution

v = b + F.*(x1 - x2);

v = regularize_mutant(v, lb, ub); % honor bounds

% binomial crossover
CR = 0.2; % tunable parameter (crossover probability) (0.2 is suggested)
crossed_over = rand(1,nvars) < CR;
new_pt = curr_point;
new_pt(crossed_over) = v(crossed_over);

end


function new_pt = polynomial_mut(curr_point, nvars, lb, ub)
%POLINOMYAL_MUT: "higly disruptive" polynomial mutation
%   some variable of the curr_point are changed according to a polynomial
%   probability distribution, to the left or to the right;
%   to have the non-higly disruptive version, delta should always be 
%   selected as min([delta1, delta2]);

Pm = 1/nvars; % probability of mutation
eta_m = 20; % tunable parameter (20 is suggested)

mutating = rand(1,nvars) < Pm;
if sum(mutating) ~= 0
    delta_1 = (curr_point(mutating) - lb(mutating))./(ub(mutating) - lb(mutating)); % normalized left variation interval
    delta_2 = (ub(mutating) - curr_point(mutating))./(ub(mutating) - lb(mutating)); % normalized right variation interval
    delta_1(ub(mutating) - lb(mutating) == 0) = 0;
    delta_2(ub(mutating) - lb(mutating) == 0) = 0;
    r = rand(1,sum(mutating));
    d1 = r <= 0.5;
    delta = delta_2;
    delta(d1) = delta_1(d1);

    delta_q = @(r, delta) (((2*r) + (1 - 2*r).*(1 - delta).^(eta_m + 1)).^(1/(eta_m + 1)) - 1).*(r <= 0.5) +...
        (1 - (2*(1 - r) + 2*(r - 0.5).*(1 - delta).^(eta_m + 1)).^(1/(eta_m + 1))).*(r > 0.5);

    dq = delta_q(r,delta);

    new_pt = curr_point;
    new_pt(mutating) = new_pt(mutating) + dq.*(ub(mutating) - lb(mutating));
    new_pt = new_pt + 0;
else
    new_pt = curr_point;
end

new_pt = regularize_mutant(new_pt, lb, ub); % honor bounds

if any(isnan(new_pt))
    disp("Hey!")
end

end


function new_pt = laplacian_mut(curr_point, nvars, lb, ub)
%LAPLACIAN_MUT: perturbs a random var by an amount sampled from Laplace
%distribution

k = randi(nvars);
b = (ub(k) - lb(k))/8; % tunable parameter
u = rand(1,1);
if u < 0.5
    incr = b*log(2*u);
else
    incr = -b*log(2 - 2*u);
end
new_pt = curr_point;
new_pt(k) = new_pt(k) + incr;

new_pt = regularize_mutant(new_pt, lb, ub); % honor bounds

end


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