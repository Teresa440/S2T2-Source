function archive = add_to_archive(archive, var, obj)
%ADD_TO_ARCHIVE: add the point specified by var and obj to the archive

id = find(archive.slot_status == 0,1,'first');
archive.vars(id, :) = var;
archive.objs(id,:) = obj;
archive.n_points = archive.n_points + 1;
archive.slot_status(id) = 1;

end