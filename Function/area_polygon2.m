function [area] = area_polygon2(P,polx,poly,polz)
nn=length(polx);
S=zeros(1,nn);
p1=P;
for i=1:1:nn-1        
    p2=[polx(i),poly(i),polz(i)];
    p3=[polx(i+1),poly(i+1),polz(i+1)];        
    S(i)=area_tri_from_point2(p1,p2,p3);
end
p3=[polx(1),poly(1),polz(1)];
p2=[polx(i+1),poly(i+1),polz(i+1)];        
S(i+1)=area_tri_from_point2(p1,p2,p3);
area=sum(S);
end