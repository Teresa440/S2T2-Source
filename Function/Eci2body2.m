function mat=Eci2body2(r_dc,angles)
phi=atan2(r_dc(2),r_dc(1))*180/pi; % azimuth [deg]
theta=asin(r_dc(3)/norm(r_dc))*180/pi; % elevation [deg]
mat=roty(90-theta)'*rotz(phi)'; % SEZ to ECI ---> v_ECI = mat*v_SEZ ### NO, hp: ECI to SEZ
mat=rotz(angles(3))*roty(angles(2))*rotx(angles(1))*mat; % ### hp: SEZ to BODY ---> v_BODY = R321*mat*v_ECI --- 321 da SEZ a BODY (?)                  
end