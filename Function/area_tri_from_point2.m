function [S] = area_tri_from_point2(p1,p2,p3)

u=[p3(1)-p1(1),p3(2)-p1(2),p3(3)-p1(3)];
v=[p2(1)-p1(1),p2(2)-p1(2),p2(3)-p1(3)];
% theta=acosd((dot(u,v))/(norm(u)*norm(v)));
% S=(norm(u)*norm(v)*sind(theta))/2;

% S=norm(cross(u,v))*0.5;

cross_uv=[u(2)*v(3)-u(3)*v(2),...
          u(3)*v(1)-u(1)*v(3),...
          u(1)*v(2)-u(2)*v(1)];
norm_uv = (cross_uv(1)^2 +...
           cross_uv(2)^2 +...
           cross_uv(3)^2)^(0.5);
S=0.5*norm_uv;

end