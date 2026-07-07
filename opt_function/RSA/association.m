function assoc_str = association(norm_arch_objs, archive_slot_status, ref_points)
%ASSOCIATION: updates line_pointer, line_dist, population_spread

assoc_str = struct;
assoc_str.line_pointer = zeros(size(norm_arch_objs,1),1); % stores which reference line is the closest to that point
assoc_str.line_dist = zeros(size(norm_arch_objs,1),1);  % store the distance from the closest reference line
assoc_str.population_spread = zeros(size(ref_points,1),1); % store how many points are associated with every reference line

for i = 1:size(norm_arch_objs,1)
    if archive_slot_status(i) ~= 0
        % distance computation is vectorized
        point = repmat(norm_arch_objs(i,:)', 1, size(ref_points,1)); % point for which distance from the reference line is computed
        line_a = zeros(size(ref_points,2), size(ref_points,1)); % first point of every reference line is the origin
        line_b = ref_points'; % second point that defines the reference line
        [min_d, id] = min(point_line_dist(point, line_a,line_b)); % find the reference line closest to "point"
        assoc_str.line_dist(i) = min_d;
        assoc_str.line_pointer(i) = id;
        assoc_str.population_spread(id) = assoc_str.population_spread(id) + 1;
    end
end

end