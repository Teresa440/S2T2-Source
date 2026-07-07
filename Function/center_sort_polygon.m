function [xyzc,xyz] = center_sort_polygon(xyz)
xyzc = mean(xyz,1);
P = xyz - xyzc;
[~,~,V] = svd(P,0);
[~,is] = sort(atan2(P*V(:,1),P*V(:,2)));
xyz = xyz(is(1:end),:);
end