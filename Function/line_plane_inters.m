function [p_end,d] = line_plane_inters(p_start,omega,n_plane,p_plane)
d=(n_plane(1)*(p_plane(1)-p_start(1))+...
   n_plane(2)*(p_plane(2)-p_start(2))+...
   n_plane(3)*(p_plane(3)-p_start(3)))/...
  (n_plane(1)*omega(1)+n_plane(2)*omega(2)+n_plane(3)*omega(3));

x=p_start(1)+omega(1)*d;
y=p_start(2)+omega(2)*d;
z=p_start(3)+omega(3)*d;

p_end=[x,y,z];
end