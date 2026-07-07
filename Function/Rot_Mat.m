function [rot] = Rot_Mat(angles)
r1=[1 0 0; 0 cosd(angles(1)) -sind(angles(1)); 0 sind(angles(1)) cosd(angles(1))];
r2=[cosd(angles(2)) 0 sind(angles(2)); 0 1 0; -sind(angles(2)) 0 cosd(angles(2))];
r3=[cosd(angles(3)) -sind(angles(3)) 0; sind(angles(3)) cosd(angles(3)) 0; 0 0 1];
rot=r3*r2*r1;
end