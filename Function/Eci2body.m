function mat=Eci2body(sat,attitude,r_dc,v_dc,cases)
r_dc_vers=r_dc/norm(r_dc);
v_dc_vers=v_dc/norm(v_dc);
phi=atan2(r_dc(2),r_dc(1))*180/pi;
theta=asin(r_dc(3)/norm(r_dc))*180/pi;

Ry = @(a) [cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)];
Rz = @(a) [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
Rx = @(a) [1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
quat2rotm = @(q) [1-2*(q(3)^2+q(4)^2), 2*(q(2)*q(3)-q(4)*q(1)), 2*(q(2)*q(4)+q(3)*q(1));
                  2*(q(2)*q(3)+q(4)*q(1)), 1-2*(q(2)^2+q(4)^2), 2*(q(3)*q(4)-q(2)*q(1));
                  2*(q(2)*q(4)-q(3)*q(1)), 2*(q(3)*q(4)+q(2)*q(1)), 1-2*(q(2)^2+q(3)^2)];

mat=Ry(90-theta)'*Rz(phi)';
switch cases
case 1
        att_info1=attitude.info1;
        att_info2=attitude.info2;
        att_nadir=attitude.nadir;
        att_vel=attitude.vel;
case 2
        att_info1=attitude.info1_cold;
        att_info2=attitude.info2_cold;
        att_nadir=attitude.nadir_cold;
        att_vel=attitude.vel_cold;
end
if strcmp(att_info1,'Nadir')
    nad=mat*(-r_dc_vers)';
    vect1=sat.geom.ext.face(att_nadir).norm;
    R=vrrotvec(vect1,nad);
    q=[cos(R(4)/2), sin(R(4)/2)*R(1:3)];
    rotm=quat2rotm(q);
    mat=rotm'*mat;
if strcmp(att_info2,'Velocity')
        vect2=sat.geom.ext.face(att_vel).norm;
        vel=mat*v_dc_vers';
        R=vrrotvec(vect2,vel);
        R(4)=R(4)*180/pi;
if isequal([1 0 0],vect1)
if R(1)<0
                R(4)=-R(4);
end
            rotm=Rx(R(4))';
elseif isequal([-1 0 0],vect1)
if R(1)>0
                R(4)=-R(4);
end
            rotm=Rx(R(4));
elseif isequal([0 1 0],vect1)
if R(2)<0
                R(4)=-R(4);
end
            rotm=Ry(R(4))';
elseif isequal([0 -1 0],vect1)
if R(2)>0
                R(4)=-R(4);
end
            rotm=Ry(R(4));
elseif isequal([0 0 1],vect1)
if R(3)<0
                R(4)=-R(4);
end
            rotm=Rz(R(4))';
elseif isequal([0 0 -1],vect1)
if R(3)>0
                R(4)=-R(4);
end
            rotm=Rz(R(4));
end
        mat=rotm*mat;
end
end
if strcmp(att_info1,'Velocity')
    vel=mat*v_dc_vers';
    vect1=sat.geom.ext.face(att_vel).norm;
    R=vrrotvec(vect1,vel);
    q=[cos(R(4)/2), sin(R(4)/2)*R(1:3)];
    rotm=quat2rotm(q);
    mat=rotm'*mat;
if strcmp(att_info2,'Nadir')
        vect2=sat.geom.ext.face(att_nadir).norm;
        nad=mat*(-r_dc_vers)';
        R=vrrotvec(vect2,nad);
        R(4)=R(4)*180/pi;
if isequal([1 0 0],vect1)
            rotm=Rx(R(4))';
elseif isequal([-1 0 0],vect1)
            rotm=Rx(R(4));
elseif isequal([0 1 0],vect1)
            rotm=Ry(R(4))';
elseif isequal([0 -1 0],vect1)
            rotm=Ry(R(4));
elseif isequal([0 0 1],vect1)
            rotm=Rz(R(4))';
elseif isequal([0 0 -1],vect1)
            rotm=Rz(R(4));
end
        mat=rotm*mat;
end
end
end