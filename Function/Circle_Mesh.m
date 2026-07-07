function [Nodes, Triangles, Quads]=Circle_Mesh(R,Nr,Nt,R_int)

if nargin < 4 || isempty(R_int)
    R_int = 0;
end

%Nodes
%------------
c=0;
for j=1:1:Nr %Number of Circles
for i=1:1:Nt %Number of Angles
c=Nt*(j-1)+i;
r_j=R_int+j*(R-R_int)/Nr;
Nodes(c,1)=r_j*cosd(360*(i-1)/Nt);
Nodes(c,2)=r_j*sind(360*(i-1)/Nt);
end
end

if R_int==0
    % Solid cylinder: single node on the axis
    Nodes(c+1,1)=0;
    Nodes(c+1,2)=0;
else
    % Hollow cylinder: inner ring of nodes at R_int
    for i=1:1:Nt
        Nodes(c+i,1)=R_int*cosd(360*(i-1)/Nt);
        Nodes(c+i,2)=R_int*sind(360*(i-1)/Nt);
    end
end

%Triangles
%------------

Triangles=[];

if R_int==0
for i=1:1:Nt-1
Triangles(i,1)=i;
Triangles(i,2)=i+1;
Triangles(i,3)=c+1;
end
Triangles(i+1,1)=i+1;
Triangles(i+1,2)=1;
Triangles(i+1,3)=c+1;
end

%Quads
%------------

Quads=[];

for j=1:1:Nr-1
for i=1:1:Nt-1
d=Nt*(j-1)+i;
Quads(d,1)=Nt*j+i;
Quads(d,2)=Nt*j+i+1;
Quads(d,3)=Nt*(j-1)+i+1;
Quads(d,4)=Nt*(j-1)+i;
end
Quads(d+1,1)=Nt*j+i+1;
Quads(d+1,2)=Nt*j+1;
Quads(d+1,3)=Nt*(j-1)+1;
Quads(d+1,4)=Nt*(j-1)+i+1;
end

if R_int>0
    % Innermost band: ring at R_int (indices c+1..c+Nt) to ring 1 (indices 1..Nt)
    Quads_in=[];
    for i=1:1:Nt-1
    Quads_in(i,1)=i;
    Quads_in(i,2)=i+1;
    Quads_in(i,3)=c+i+1;
    Quads_in(i,4)=c+i;
    end
    Quads_in(i+1,1)=i+1;
    Quads_in(i+1,2)=1;
    Quads_in(i+1,3)=c+1;
    Quads_in(i+1,4)=c+i+1;

    Quads=[Quads_in; Quads];
end

end