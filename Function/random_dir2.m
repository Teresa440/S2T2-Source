function [omega] = random_dir2(R,omega)

theta=asin((rand(1))^(1/2));
phi=2*pi*rand(1);

% theta = acos(rand(1));
% phi = 2*pi*rand(1);

sin_theta = sin(theta);
omega(1) = sin_theta*cos(phi);
omega(2) = sin_theta*sin(phi);
omega(3) = cos(theta);

omega=R*omega;
end