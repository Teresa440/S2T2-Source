% Following functions imported from:

% ----------------------------------------------------------------------- %
% Example of use of the funcion MOPSO.m, which performs a Multi-Objective %
% Particle Swarm Optimization (MOPSO), based on Coello2004.               %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez Cagigal                                      %
%   Date:    15/03/2017                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%       Coello, C. A. C., Pulido, G. T., & Lechuga, M. S. (2004). Handling%
%       multiple objectives with particle swarm optimization. IEEE Tran-  %
%       sactions on evolutionary computation, 8(3), 256-279.              %
% ----------------------------------------------------------------------- %

function d = dominates(x,y)
%DOMINATES: returns 1 if x dominates y and 0 otherwise

d = all(x<=y,2) & any(x<y,2);

end