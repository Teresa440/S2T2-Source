function [r_out,v_out,T_p]=orb_calc(env,orbit,sim_data)
%la funzione ha come output:
% -r= coordinate orbita in sistema di riferimento Geocentrico per n_orb
% orbite
% -T_p=Periodo orbitale di un orbita
if nargin<3
    warning('Insert the following structure: environment, orbit, sim_data')
else
    Rp=env.Rp;
    jdate=sim_data.jdat;
    [~, decl, ~] = sun2 (jdate);
    %n_orb=sim_data.n_orbit;
    step=sim_data.step;
    if strcmp(orbit.par_orb,'Basic')==1
        beta_angle=orbit.beta_angle;
        z=orbit.z;
        %L1=[1 0 0;0 cos(eps) -sin(eps);0 sin(eps) cos(eps)];
        incl=beta_angle+decl*180/pi;        
        e=0;
        a=Rp+z;                                                            %[km]
        omega=0;
        OMEGA=0;
        
    elseif strcmp(orbit.par_orb,'Classic')==1
        e=orbit.eccentricity;
        a=orbit.SA;
        omega=orbit.AP;                                                    %argument Periapsid
        OMEGA=orbit.LAN;
        incl=orbit.inclination;
    end
    
    %tracciare l'orbita
    mu=env.mu;
    T_p=round(2*pi*(a^3/mu)^0.5);
    tempo(T_p);
    %T=linspace(0,n_orb*T_p,(n_orb*T_p)/step);
    T=linspace(0,T_p,T_p/step);
    r=zeros(length(T),3);
    E=zeros(length(T),1);
    ni_polar=zeros(length(T),1);
    r_polar=zeros(length(T),1);
    r_perifocale=zeros(length(T),1);
    for i=1:length(T)
        [E(i,1),ni_polar(i,1)]=inverseproblem(T(i),e,a,mu);                %ni polar esce in deg
        r_polar(i,1)=(a*(1-e^2))/(1+e*cos(ni_polar(i,1)*pi/180));
        % vettore r nel piano perifocale:
        r_perifocale(i,1)=r_polar(i,1)*cos(ni_polar(i,1)*pi/180);
        r_perifocale(i,2)=r_polar(i,1)*sin(ni_polar(i,1)*pi/180);
        r_perifocale(i,3)=0;
        v_perifocale(i,1)=(mu/(a*(1-e^2)))^1/2*(-sin(ni_polar(i,1)*pi/180));
        v_perifocale(i,2)=(mu/(a*(1-e^2)))^1/2*(e+cos(ni_polar(i,1)*pi/180));
        v_perifocale(i,3)=0;
        r(i,:)=per_ijk(r_perifocale(i,:)',OMEGA,omega,incl);
        v(i,:)=per_ijk(v_perifocale(i,:)',OMEGA,omega,incl);
    end
    r_out=[r T'];
    v_out=[v T'];
end
end