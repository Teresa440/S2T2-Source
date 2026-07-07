% Following functions imported from:

% ----------------------------------------------------------------------- %
% Example of use of the funcion MOPSO.m, which performs a Multi-Objective %
% Particle Swarm Optimization (MOPSO), based on Coello2004.               %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez Cagigal                                      %
%   Date:    15/03/2017                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%       Coello, C. A. C., Pulido, G. T., & Lechuga, M. S. (2004). Handling%
%       multiple objectives with particle swarm optimization. IEEE Tran-  %
%       sactions on evolutionary computation, 8(3), 256-279.              %
% ----------------------------------------------------------------------- %

function dom_vector = check_domination(fitness)
%CHECK_DOMINATION: checks the domination between the population, it
% returns a vector that indicates if each particle is dominated (1) or not

Np = size(fitness,1);
dom_vector = zeros(Np,1);
if Np > 1 % added this to consider the case Np = 1
    all_perm = nchoosek(1:Np,2); % possible permutations
    all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];
    d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
    dominated_particles = unique(all_perm(d==1,2));
    dom_vector(dominated_particles) = 1;
end

end

