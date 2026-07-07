function [r]=per_ijk(R,OMEGA,omega,incl)
    L1=[cos(OMEGA*pi/180) -sin(OMEGA*pi/180) 0;
    sin(OMEGA*pi/180) cos(OMEGA*pi/180) 0;
    0 0 1];

L2=[1 0 0;
    0 cos(incl*pi/180) -sin(incl*pi/180);
    0 sin(incl*pi/180) cos(incl*pi/180)];

L3=[cos(omega*pi/180) -sin(omega*pi/180) 0;
    sin(omega*pi/180) cos(omega*pi/180) 0;
    0 0 1];
L=L1*L2*L3;
r=L*R;
return