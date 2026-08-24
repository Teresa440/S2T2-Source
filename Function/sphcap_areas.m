function [Afb,Aft,Afl,Acb,Act,Acl,V,Geom_r] = sphcap_areas(R_int,R_out,h,Ntheta,Nr_cap,Nt)


dphi=2*pi/Nt;
theta_max=pi/2;

Rs_bound=zeros(1,Nr_cap+1);
for k=0:1:Nr_cap
    Rs_bound(k+1)=R_int+k*(R_out-R_int)/Nr_cap; % Rs=r_rim 
end

Afb=zeros(Ntheta,Nr_cap,6);
Aft=zeros(Ntheta,Nr_cap,6);
Afl=zeros(Ntheta,Nr_cap,6);
Acb=zeros(Ntheta,Nr_cap,6);
Act=zeros(Ntheta,Nr_cap,6);
Acl=zeros(Ntheta,Nr_cap,6);
V=zeros(Ntheta,Nr_cap);
Geom_r=zeros(Ntheta,Nr_cap);

for hh=1:1:Nr_cap

r_i=Rs_bound(hh);
r_o=Rs_bound(hh+1);
t=r_o-r_i;

    for j=1:1:Ntheta

    theta_lo=(j-1)*theta_max/Ntheta;
    theta_hi=j*theta_max/Ntheta;
    dtheta=theta_max/Ntheta;

    A(1)=0.5*sin(dtheta)*(r_o^2-r_i^2);                     % +phi
    A(4)=A(1);                                              % -phi
    cos_gamma_hi=sin(theta_hi)^2*cos(dphi)+cos(theta_hi)^2;
    cos_gamma_lo=sin(theta_lo)^2*cos(dphi)+cos(theta_lo)^2;
    A(2)=0.5*sin(acos(cos_gamma_hi))*(r_o^2-r_i^2);         % +theta 
    A(5)=0.5*sin(acos(cos_gamma_lo))*(r_o^2-r_i^2);         % -theta 

    A(3)=quad_area_fan(pos3(r_o,theta_lo,0),pos3(r_o,theta_lo,dphi),pos3(r_o,theta_hi,dphi),pos3(r_o,theta_hi,0));
    A(6)=quad_area_fan(pos3(r_i,theta_lo,0),pos3(r_i,theta_lo,dphi),pos3(r_i,theta_hi,dphi),pos3(r_i,theta_hi,0));

    V(j,hh)=(r_o^3-r_i^3)/3*(cos(theta_lo)-cos(theta_hi))*dphi;

    dOmega=dphi*(cos(theta_lo)-cos(theta_hi));
    Geom_r(j,hh)=dOmega*r_i*r_o/t;

    Afb_row=zeros(1,6);
    Aft_row=zeros(1,6);
    Afl_row=zeros(1,6);
    Acb_row=A;
    Act_row=A;
    Acl_row=A;

    if hh==1
        Afb_row(6)=A(6);
        Acb_row(6)=0;
    end
    if hh==Nr_cap
        Aft_row(3)=A(3);
        Act_row(3)=0;
    end
    if j==Ntheta
        Afb_row(2)=A(2);
        Aft_row(2)=A(2);
        Afl_row(2)=A(2);
        Acb_row(2)=0;
        Act_row(2)=0;
        Acl_row(2)=0;
    end

    Afb(j,hh,:)=Afb_row;
    Aft(j,hh,:)=Aft_row;
    Afl(j,hh,:)=Afl_row;
    Acb(j,hh,:)=Acb_row;
    Act(j,hh,:)=Act_row;
    Acl(j,hh,:)=Acl_row;

    end
end

end

function P = pos3(r,theta,phi)
P = r*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
end

function A = quad_area_fan(P1,P2,P3,P4)

P=(P1+P2+P3+P4)/4;
A = 0.5*norm(cross(P2-P,P1-P)) + 0.5*norm(cross(P3-P,P2-P)) ...
  + 0.5*norm(cross(P4-P,P3-P)) + 0.5*norm(cross(P1-P,P4-P));
end
