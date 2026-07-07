function [R] = rot_mat_emiss(n,c,p)
t1=(p-c)/norm(p-c);
t2=cross(n,t1);

A=[t1',t2',n'];

R=A;


end