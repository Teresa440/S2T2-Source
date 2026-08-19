function [Afb,Aft,Afl,Acb,Act,Acl,V,Geom_r] = sphcap_areas(R_int,R_out,h,Ntheta,Nr_cap,Nt)

Rs_out=(R_out^2+h^2)/(2*h);
dphi=2*pi/Nt;

Rs_bound=zeros(1,Nr_cap+1);
theta_max_bound=zeros(1,Nr_cap+1);
for k=0:1:Nr_cap
    r_rim=R_int+k*(R_out-R_int)/Nr_cap;
    Rs_bound(k+1)=sqrt(r_rim^2+(Rs_out-h)^2);
    theta_max_bound(k+1)=atan2(r_rim,Rs_out-h);
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

theta_max_i=theta_max_bound(hh);
theta_max_o=theta_max_bound(hh+1);

% theta_max varia col layer radiale (bordo piatto, vedi Sph_Shell_Mesh.m):
% la cella (j,hh) non e' un elemento sferico puro. Usato solo per A_phi/
% A_theta/V (le facce radiali A(3)/A(6) restano esatte, ciascuna sulla
% propria sfera theta_max_o/theta_max_i) e per Geom_r (vedi sotto). E'
% un errore di discretizzazione: deve annullarsi aumentando Nr_cap
% (verifica di convergenza in Fase 9).
theta_max_avg=0.5*(theta_max_i+theta_max_o);

    for j=1:1:Ntheta

    theta_lo_i=(j-1)*theta_max_i/Ntheta;
    theta_hi_i=j*theta_max_i/Ntheta;
    theta_lo_o=(j-1)*theta_max_o/Ntheta;
    theta_hi_o=j*theta_max_o/Ntheta;
    theta_lo_avg=(j-1)*theta_max_avg/Ntheta;
    theta_hi_avg=j*theta_max_avg/Ntheta;
    dtheta_avg=theta_max_avg/Ntheta;

    % Aree Af/Ac: pannello PIATTO vero (poligono reale della mesh), non
    % superficie curva ideale -- stessa convenzione di cylinder_areas.m
    % (verificato: la sua A(2) coincide esattamente con l'area del quad
    % 3D reale della mesh, non con l'arco r*dtheta; differenza ~1% su
    % una mesh di riferimento). Af/Ac alimentano irraggiamento (Af_tot)
    % ed eventuali direzioni di conduzione generiche (theta/phi): DEVONO
    % essere il pannello vero, non la fisica a sfera liscia.
    %
    % Facce phi (1,4) e theta (2,5): dimostrabilmente PIANE (ciascuna
    % giace su due rette per il centro delle sfere -> definiscono un
    % piano), quindi hanno una formula chiusa esatta -- verificata contro
    % la triangolazione diretta dei vertici (differenza ~1e-13).
    A(1)=0.5*sin(dtheta_avg)*(r_o^2-r_i^2);                 % +phi
    A(4)=A(1);                                              % -phi
    cos_gamma_hi=sin(theta_hi_avg)^2*cos(dphi)+cos(theta_hi_avg)^2;
    cos_gamma_lo=sin(theta_lo_avg)^2*cos(dphi)+cos(theta_lo_avg)^2;
    A(2)=0.5*sin(acos(cos_gamma_hi))*(r_o^2-r_i^2);         % +theta (verso il bordo)
    A(5)=0.5*sin(acos(cos_gamma_lo))*(r_o^2-r_i^2);         % -theta (verso l'apice; =0 esatto per j=1)

    % Facce radiali (3,6): NON piane in generale (i 4 vertici stanno
    % sulla stessa sfera, non su rette per il centro) -- serve la
    % triangolazione a ventaglio dal centroide, stesso schema gia' usato
    % da area_polygon2.m/area_tri_from_point2.m per pannelli non piani.
    % Theta esatto del proprio layer (o/i), non mediato -- confermato
    % corretto in precedenza, nessun compromesso qui.
    A(3)=quad_area_fan(pos3(r_o,theta_lo_o,0),pos3(r_o,theta_lo_o,dphi),pos3(r_o,theta_hi_o,dphi),pos3(r_o,theta_hi_o,0));
    A(6)=quad_area_fan(pos3(r_i,theta_lo_i,0),pos3(r_i,theta_lo_i,dphi),pos3(r_i,theta_hi_i,dphi),pos3(r_i,theta_hi_i,0));

    V(j,hh)=(r_o^3-r_i^3)/3*(cos(theta_lo_avg)-cos(theta_hi_avg))*dphi;

    % Fattore geometrico per la conduttanza radiale in forma chiusa
    % (Fase 5, TMM2.m): G_r_cella = k*dOmega_avg*r_i*r_o/(r_o-r_i), stesso
    % theta_max medio di A_theta/A_phi/V (non quello esatto di A(3)/A(6)).
    % Convergenza verificata separatamente con un integrale 1D indipendente,
    % sia sul caso aggregato (apice) sia su una cella theta generica:
    % errore -> 0 come ~1/Nr_cap^2 (secondo ordine, pulito) in entrambi i casi.
    dOmega_avg=dphi*(cos(theta_lo_avg)-cos(theta_hi_avg));
    Geom_r(j,hh)=dOmega_avg*r_i*r_o/t;

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
% Area di un quadrilatero (eventualmente non piano) per triangolazione a
% ventaglio dal centroide -- stesso schema di area_polygon2.m.
P=(P1+P2+P3+P4)/4;
A = 0.5*norm(cross(P2-P,P1-P)) + 0.5*norm(cross(P3-P,P2-P)) ...
  + 0.5*norm(cross(P4-P,P3-P)) + 0.5*norm(cross(P1-P,P4-P));
end
