function d = point_line_dist(point, line_a, line_b)
%POINT_LINE_DIST: distance from point to line, works in 
%N-dimensions
%   inputs are supposed to be column vectors, or many column vectors next
%   to each other, for vectorized computation

t = (dot(point, line_b - line_a) - dot(line_a, line_b - line_a))./(vecnorm(line_b - line_a).^2);
d = vecnorm(point - (line_a + t.*(line_b - line_a)));

end