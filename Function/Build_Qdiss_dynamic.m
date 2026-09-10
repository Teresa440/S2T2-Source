function [Qdiss_dyn] = Build_Qdiss_dynamic(Q_diss_modes, mode_schedule, dt, r_dc, sim_data)

%warning: each orbit is assumed to have the same duration and the same mode schedule. If this is not the case, the code needs to be modified

n_nodes = size(Q_diss_modes, 1);
n_per_orbit = length(r_dc) - 1;

% only t_start matters: the active mode is the last one that has started
% by time t, so t_end is purely informational and gaps/mismatches between
% a mode's t_end and the next mode's t_start can't break the schedule
[t_starts, order] = sort(mode_schedule(:,1));
Q_diss_modes_sorted = Q_diss_modes(:, order);

Qdiss_orbit = zeros(n_nodes, n_per_orbit);

for k = 1:n_per_orbit
    t = (k-1)*dt;   % time within the orbit, orbit starts at t=0
    idx_mode = find(t_starts <= t, 1, 'last');

    if isempty(idx_mode)
        error('No operational mode starts at or before t = %g s (first mode must start at t=0)', t);
    end

    Qdiss_orbit(:, k) = Q_diss_modes_sorted(:, idx_mode);
end

Qdiss_dyn = repmat(Qdiss_orbit, 1, sim_data.n_orbit);

end