function [archive] = initialize_archive(fun, options, nvars, nobjs, lb, ub, ref_points)
%INITIALIZE_ARCHIVE: create the initial archive for RSA
%   gamma*SL random points are generated, then only non-dominated point are
%   kept.
%   Dominance problem could be sped up using "ND-Tree-based update: a Fast 
%   Algorithm for the Dynamic Non-Dominance Problem"
%   (https://arxiv.org/pdf/1603.04798.pdf)

max_iter = 10; % maximum number of iteration for each point to refine (optional) the initial solutions

archive = struct;
% random initial population
archive.vars = lb + (ub - lb).*rand(options.gamma*options.SL,nvars); % decision variables corresponding to points of the archive
archive.objs = zeros(options.gamma*options.SL,nobjs); % objectives values corresponding to points of the archive
archive.slot_status = zeros(options.gamma*options.SL,1); % archive.slot_status contains 1 if it is occupied, 0 if it is not, 2 if it is non-dominated
archive.pt_id = (1:options.gamma*options.SL)'; % index of the points in the archive

for i = 1:size(archive.vars,1)
    pt_val = fun(archive.vars(i,:)); % objectives values corresponding to points of the archive
    if size(pt_val,1) > size(pt_val,2)
        pt_val = pt_val'; % if fun outputs a col vector it is transposed to a row vector
    end
    archive.objs(i,:) = pt_val;
end

% Local optimization of every variable (optional)
if options.do_initial_local_search == true
    archive = local_optimization(options, archive, fun, max_iter, ub, lb, nvars);
end

dominated_points = check_domination(archive.objs); % lists the dominated points solving the dominance problem
archive.slot_status = dominated_points; % dominated points
archive.slot_status(archive.slot_status == 0) = 2; % non-dominated points
archive.n_points = sum(archive.slot_status ~= 0); % number of points, archive size

% only non-dominated points are kept at the end of initialization
archive = remove_from_archive(archive.slot_status == 1,archive);

if archive.n_points > options.HL
    % if front points exceed HL (improbable) they are clusterd down to HL
    archive = cluster(archive, ref_points, options, nobjs);
end

end

%% Local functions

function archive = local_optimization(options, archive, fun, max_iter, ub, lb, nvars)
%LOCAL_OPTIMIZATION: tries to improve initial point by using a local search

for i = 1:size(archive.vars,1)
    curr_var = archive.vars(i,:);
    curr_obj = fun(curr_var);
    for k = 1:max_iter
        new_pt = perturb(archive, options, nvars, curr_var, 1, lb, ub); % perturb the current point
        new_obj = fun(new_pt);
        if size(new_obj,1) > size(new_obj,2)
            new_pt_val = new_pt_val'; % if fun outputs a col vector it is transposed to a row vector
        end
        if dominates(new_obj, curr_obj) % new point is kept only if it dominatesn the old one
            curr_var = new_var;
            curr_obj = new_obj;
        end
    end
    archive.vars(i,:) = curr_var;
    archive.objs(i,:) = curr_obj;
end

end