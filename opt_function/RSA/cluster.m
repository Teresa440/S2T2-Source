function archive = cluster(archive, ref_points, options, nobjs)
%CLUSTER Cluster points to reduce archive size back to HL (hard limit)
%   adaptive normalization is performed first, then the most distant point
%   from the line with higher population_spread are removed

norm_arch_objs = adaptive_normalization(archive,nobjs);
assoc_str = association(norm_arch_objs, archive.slot_status, ref_points);

% clustering
while archive.n_points > options.HL
    [~, ref_line_id] = max(assoc_str.population_spread); % reference line with maximum population spread
    point_list = (assoc_str.line_pointer == ref_line_id); % logical array of points associated to ref_line_id
    [~, pt_number] = max(assoc_str.line_dist(point_list));
    temp_pt_id = archive.pt_id(point_list);
    delete_id = temp_pt_id(pt_number); % id of the point to be deleted from the archive
    % deleting the point:
    archive = remove_from_archive(delete_id,archive);
    assoc_str.line_pointer(delete_id) = 0;
    assoc_str.line_dist(delete_id) = 0;
    assoc_str.population_spread(ref_line_id) = assoc_str.population_spread(ref_line_id) - 1;
end

end

%% Local functions

function plot_normalized_archive(nobjs, norm_arch_objs, ref_points)
%PLOT_NORMALIZED_ARCHIVE: debug function used to see the clustering process
%in 2 or 3 dimensions

if nobjs == 2
    figure
    plot(norm_arch_objs(norm_arch_objs(:,1)~=0,1), norm_arch_objs(norm_arch_objs(:,2)~=0,2),"pentagram")
    grid on
    hold on
    for i = 1:size(ref_points,1)
        plot([0 ref_points(i,1)], [0 ref_points(i,2)],'k')
        plot([ref_points(i,1)], [ref_points(i,2)],'ok')
    end
elseif nobjs == 3
    figure
    plot3(norm_arch_objs(norm_arch_objs(:,1)~=0,1), norm_arch_objs(norm_arch_objs(:,2)~=0,2),norm_arch_objs(norm_arch_objs(:,3)~=0,3),"pentagram")
    grid on
    hold on
    for i = 1:size(ref_points,1)
        plot3([0 ref_points(i,1)], [0 ref_points(i,2)], [0 ref_points(i,3)],'k')
        plot3([ref_points(i,1)], [ref_points(i,2)], [ref_points(i,3)],'ok')
    end
end

end