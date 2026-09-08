function [Qdiss_dyn] = Build_Qdiss_dynamic(Q_diss_modes, mode_schedule, dt, r_dc, sim_data)

n_nodes = size(Q_diss_modes, 1);
N_modes = size(mode_schedule, 1);

% instants per orbit: same exact value used everywhere else (select_trans3,
% transient, transient3_pruning) to size T/Tempo, so it stays consistent
% with the rest of the pipeline without depending on Tempo being filled in yet
n_per_orbit = length(r_dc) - 1;

Qdiss_orbit = zeros(n_nodes, n_per_orbit);

for k = 1:n_per_orbit
    t = (k-1)*dt;   % time within the orbit, orbit starts at t=0
    idx_mode = [];

    for m = 1:N_modes
        is_last_mode = (m == N_modes);
        if is_last_mode
            in_range = (t >= mode_schedule(m,1)) && (t <= mode_schedule(m,2));
        else
            in_range = (t >= mode_schedule(m,1)) && (t < mode_schedule(m,2));
        end
        if in_range
            idx_mode = m;
            break
        end
    end

    if isempty(idx_mode)
        error('There is no operational mode defined for t = %g s', t);
    end

    Qdiss_orbit(:, k) = Q_diss_modes(:, idx_mode);
end

Qdiss_dyn = repmat(Qdiss_orbit, 1, sim_data.n_orbit);

end