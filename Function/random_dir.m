function [omega] = random_dir(R)

theta=asin((rand(1))^(1/2));
phi=2*pi*rand(1);

% theta = acos(rand(1));
% phi = 2*pi*rand(1);                

omegaG=[sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];

omega=R*omegaG;
end