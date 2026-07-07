function [Afb,Aft,Afl,Acb,Act,Acl,V] = cylinder_areas(R,Nr,Nt,L,Nz,R_int)

if nargin < 6 || isempty(R_int)
    R_int = 0;
end

a=(R-R_int)/Nr;
alfa=(360/Nt)/2;
dz=L/(Nz-1);

Afb=zeros(Nr,6);
Aft=zeros(Nr,6);
Afl=zeros(Nr,6);
Acb=zeros(Nr,6);
Act=zeros(Nr,6);
Acl=zeros(Nr,6);

V=zeros(1,Nr);


for j=1:1:Nr
    r_out=R_int+j*a;
    r_in=R_int+(j-1)*a;

    A(1)=a*dz;
    A(2)=2*r_out*sind(alfa)*dz;
    A(3)=0.5*sind(2*alfa)*(r_out^2-r_in^2);
    A(4)=a*dz;
    A(5)=2*r_in*sind(alfa)*dz;
    A(6)=A(3);
    V(j)=A(3)*dz;

    Afb(j,:)=[0 0 0 0 0 A(6)];
    Aft(j,:)=[0 0 A(3) 0 0 0];
    Afl(j,:)=[0 0 0 0 0 0];

    Acb(j,:)=[A(1) A(2) A(3) A(4) A(5) 0];
    Act(j,:)=[A(1) A(2) 0 A(4) A(5) A(6)];
    Acl(j,:)=[A(1) A(2) A(3) A(4) A(5) A(6)];

end

% Outer lateral surface (ring Nr): the +r face is the real external
% surface of the cylinder, not a conduction contact.
Afb(Nr,2)=A(2);
Aft(Nr,2)=A(2);
Afl(Nr,2)=A(2);

Acb(Nr,2)=0;
Act(Nr,2)=0;
Acl(Nr,2)=0;

if R_int>0
    % Inner bore surface (ring 1): the -r face is the real internal
    % surface of the hollow cylinder, not a conduction contact towards
    % a non-existent ring 0.
    a_in=2*R_int*sind(alfa)*dz;

    Afb(1,5)=a_in;
    Aft(1,5)=a_in;
    Afl(1,5)=a_in;

    Acb(1,5)=0;
    Act(1,5)=0;
    Acl(1,5)=0;
end


