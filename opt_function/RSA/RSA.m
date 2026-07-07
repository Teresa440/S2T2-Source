%___________________________________________________________________%
% Reference Point Based Archived Many Objective Simulated Annealing %
%                                                                   %
%  Developed in MATLAB R2022b                                       %
%                                                                   %
%  Author and programmer: Davide Cosenza                            %
%                                                                   %
%         e-Mail: davidecosenza98@gmail.com                         %
%                 davide.cosenza@studenti.polito.it                 %
%                                                                   %
%                                                                   %
%    Code adaptation based on the main paper:                       %
%                                                                   %
%        Raunak Sengupta, Sriparna Saha,                            %
%        Reference point based archived                             %
%        many objective simulated annealing,                        %
%        Information Sciences, Volume 467, 2018, Pages 725-749,     %
%        ISSN 0020-0255,                                            %
%         https://doi.org/10.1016/j.ins.2018.05.013.                %
%                                                                   %
%___________________________________________________________________%

% I acknowledge that this version of RSA has been written using
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_pareto, fval_pareto, x, fval, archive, options] = RSA(fun, lb, ub, options)

%RSA: Reference Point Based Archived Many Objective Simulated Annealing
% 
%   Finds Pareto front of solutions which minimizes "fun" using RSA 
%   algorithm. 
%   "fun" can be a n-inputs m-outputs function, with n = number of
%   variables ("nvars" in the code) and m = number of objective ("nobjs" in
%   the code).
%   Function inputs must be bounded by lower bounds specified by a 
%   row vector "lb" and upper bounds specified by a row vector "ub".
%   Returns the variable values on the Pareto front "x_pareto" with 
%   relative fitnesses "fval_pareto", then returns the complete collection
%   of variable and fitnesses "x" and "fval" and finally returns the
%   archive at the end of the optimization and the options used,
%   "archive" and "options", respectively.
%
%   Code could be further improved by using vectorization of fun evaluation

if size(lb,1) ~= 1 || size(ub,1) ~= 1 % check dimensions of lb and ub
    error("Error: lb and ub must be row vectors")
end

nvars = length(lb); % number of variables
nobjs = length(fun(lb)); % number of objectives

if ~exist('options','var')
    options = struct;
    options.gamma = 2; % multiplication factor of initial archive size
    options.HL = 100; % hard limit of the archive
    options.SL = 200; % soft limit of the archive (suggested: 2 or 3 times HL) 
    options.Tmax = 100; % initial tempereature
    options.Tmin = 0.000001; % final tempereature
    options.iter = 100; % iterations at each temperature
    options.alpha = 0.99; % factor of decrease in temperature, cooling rate (suggested: 0.99) (AMOSA uses 0.8) % 0.9908319449
    options.max_fe = 10000;
    options.do_initial_local_search = false; % (suggested: false)
    if nobjs >= 2 && nobjs <= 15
        suggested_ref_points = [0 100 12 8 6 4 4 3 3 3 2 2 2 2 2];
        options.ref_point_per_dimension = suggested_ref_points(nobjs); % number of reference point in each axis of the nobjs dimensions
    end
    % perturbation options
    options.flag1 = 1; % if set to 1 polynomial mutation is performed after differential mutation
    options.prob1 = 0.1; % probability of performin laplacian mutation after differential mutation
    options.flag2 = 1; % if set to 1 polynomial mutation is performed after Simulated Binary Crossover (SBX) mutation
    options.prob2 = 0.1; % probability of performin laplacian mutation after Simulated Binary Crossover (SBX) mutation
    options.switch_factor = 0.5; % fraction of iteration where Simulated Binary Crossover (SBX) mutation is performed, in the other iteration differential mutation is used instead
    % behaviour options
    options.plot_results = false;
    options.plot_rejected = false;
    options.verbose = true;
end

if options.verbose
    total_fun_evaluation = round(log(options.Tmin/options.Tmax)/log(options.alpha))*options.iter;
    total_fun_evaluation = total_fun_evaluation + options.SL*options.gamma;
    if options.do_initial_local_search
        total_fun_evaluation = total_fun_evaluation + 10*options.SL*options.gamma;
    end
    fprintf("Total number of function evaluation expected: %d\n\n", total_fun_evaluation);
end

ref_points = gen_ref_points(options.ref_point_per_dimension, nobjs); % generates the reference point (lines)

archive = initialize_archive(fun, options, nvars, nobjs, lb, ub, ref_points); % creates the starting archive
f_e = options.SL*options.gamma;

T = options.Tmax;
cont = 0; % counts up every time T changes
rejected = 0; % number of times new_pt gets rejected
rejected_iter = 0;

while T > options.Tmin && f_e < options.max_fe
    for i = 1:options.iter
        pt_list = archive.slot_status ~= 0;
        arch_id = archive.pt_id(pt_list);
        id_c = randi(archive.n_points); % number of point selected
        a_id_c = arch_id(id_c); % id of the row of the archive corresponding to the selected point
        new_pt = perturb(archive, options, nvars, archive.vars(a_id_c,:), i, lb, ub); % point is perturbed
        new_pt_val = fun(new_pt);
        f_e = f_e + 1;
        if size(new_pt_val,1) > size(new_pt_val,2)
            new_pt_val = new_pt_val'; % if fun outputs a col vector it is transposed to a row vector
        end
        new_pt_val_rep = repmat(new_pt_val,archive.n_points, 1);
        l = dominates(new_pt_val_rep, archive.objs(pt_list,:)); % number of times new_pt dominates the points in the archive
        k = dominates(archive.objs(pt_list,:), new_pt_val_rep); % number of times the points in the archive dominate new_pt
        prob = exp((sum(l) - sum(k))/(archive.n_points*T)); % acceptance probability
        if dominates(archive.objs(a_id_c,:), new_pt_val) % Case 1: the current point dominates the new point
            % the new point is added to the archive with a probability prob
            if rand(1,1) < prob
                archive = add_to_archive(archive, new_pt, new_pt_val);
                if archive.n_points > options.SL
                    archive = cluster(archive, ref_points, options, nobjs);
                end
            end
        elseif sum(l) > sum(k) % Case 2a: the new point dominates the current point with sum(l) > sum(k)
            % the new point is added to the archive and the points in the archive dominated by it are removed
            archive = add_to_archive(archive, new_pt, new_pt_val);
            archive = remove_from_archive(arch_id(l), archive);
            if archive.n_points > options.SL
                archive = cluster(archive, ref_points, options, nobjs);
            end
        elseif rand(1,1) < prob % Case 2b: the new point dominates the current point with sum(l) <= sum(k)
            % the new point is added to the archive with a probability prob
            archive = add_to_archive(archive, new_pt, new_pt_val);
            if archive.n_points > options.SL
                archive = cluster(archive, ref_points, options, nobjs);
            end
        else
            % counting the rejected points
            rejected = rejected + 1;
            rejected_iter = rejected_iter + 1;
        end
    end

    T = T*options.alpha; % temperature update
    cont = cont + 1; % number of T changes

    % plots and verbose
    if options.plot_results || options.verbose || options.plot_rejected
        cont_skip = 50;
        if rem(cont,cont_skip) == 0 && options.plot_rejected
            figure(6)
            plot(cont,rejected_iter/cont_skip,'.k')
            rejected_iter = 0;
            grid on
            hold on
            title("Moving average of the number rejected point")
            xlabel("Temperature changes")
            ylabel("Number of rejected points")
        end
        if rem(cont,cont_skip) == 0 && options.plot_results
            archive = verbose_plot(archive, options, cont, T, ref_points);
        end
    end

end

archive = cluster(archive, ref_points, options, nobjs); % final clustering

% final plots and verbose
if options.plot_results || options.verbose
    archive = verbose_plot(archive, options, cont, T, ref_points);
end

pt_list = archive.slot_status ~= 0;
dominated_points = check_domination(archive.objs(pt_list,:)); % determines domination
dominated_points(dominated_points == 0) = 2;
archive.slot_status(pt_list) = dominated_points;
x_pareto = archive.vars(archive.slot_status == 2,:); % variables of the points on the front
fval_pareto = archive.objs(archive.slot_status == 2,:); % fitnesses of the points on the front
x = archive.vars(pt_list,:); % return value: variables values of the points in the archive
fval = archive.objs(pt_list,:); % return value: variables fitness value of the points in the archive

end

%% Local functions

function archive = verbose_plot(archive, options, cont, T, ref_points)
if options.plot_results
    pt_list = archive.slot_status ~= 0;
    dominated_points = check_domination(archive.objs(pt_list,:));
    dominated_points(dominated_points == 0) = 2;
    archive.slot_status(pt_list) = dominated_points;
    plot_archive(archive.objs, archive.slot_status, ref_points, options.ax,  cont, T, options.Tmax, options.Tmin)
    % axis([0 5 0 5 0 5])
    % axis equal
end
if options.verbose
    fprintf("Cont: %d,\tT = %f, Archive size = %d\n",cont,T,archive.n_points);
end

end