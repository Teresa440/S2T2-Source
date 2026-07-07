function [Afb,Aft,Afl,Acb,Act,Acl,V] = cylinder_areas(R,Nr,Nt,L,Nz)

a=R/Nr;
alfa=(360/Nt)/2;
c=a*cosd(alfa);
b=a*sind(alfa);
dz=L/(Nz-1);

Afb=zeros(Nr,6);
Aft=zeros(Nr,6);
Afl=zeros(Nr,6);
Acb=zeros(Nr,6);
Act=zeros(Nr,6);
Acl=zeros(Nr,6);

V=zeros(1,Nr);


for j=1:1:Nr
    A(1)=a*dz;
    A(2)=(2*b*dz)*j;
    A(3)=(c*b)+2*(c*b)*(j-1);
    A(4)=a*dz;
    A(5)=(2*b*dz)*(j-1);
    A(6)=(c*b)+2*(c*b)*(j-1);
    V(j)=A(3)*dz; 

    Afb(j,:)=[0 0 0 0 0 A(6)];
    Aft(j,:)=[0 0 A(3) 0 0 0];
    Afl(j,:)=[0 0 0 0 0 0];

    Acb(j,:)=[A(1) A(2) A(3) A(4) A(5) 0];
    Act(j,:)=[A(1) A(2) 0 A(4) A(5) A(6)];
    Acl(j,:)=[A(1) A(2) A(3) A(4) A(5) A(6)];

end

Afb(j,2)=A(2);
Aft(j,2)=A(2);
Afl(j,2)=A(2);

Acb(j,2)=0;
Act(j,2)=0;
Acl(j,2)=0;


