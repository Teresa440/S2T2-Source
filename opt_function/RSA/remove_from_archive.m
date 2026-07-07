function archive = remove_from_archive(rows,archive)
%REMOVE_FROM_ARCHIVE: removes rows specified by rows in every variable of
%the archive, and substitute them with zeros.
%   rows must be a logical array or an array of integers

archive.vars(rows,:) = 0;
archive.objs(rows,:) = 0;
archive.slot_status(rows) = 0;
% archive.n_points = sum(archive.slot_status ~= 0);
if islogical(rows)
    archive.n_points = archive.n_points - sum(rows);
else
    archive.n_points = archive.n_points - length(rows);
end

end

