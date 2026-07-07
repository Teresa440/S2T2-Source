function [normal] = normal_from_points(xyz)
center=[mean(xyz(:,1)),mean(xyz(:,2)),mean(xyz(:,3))];

u=[xyz(1,1)-center(1),xyz(1,2)-center(2),xyz(1,3)-center(3)];
v=[xyz(2,1)-center(1),xyz(2,2)-center(2),xyz(2,3)-center(3)];

normal_vect=cross(u,v);
normal=normal_vect/norm(normal_vect);
end